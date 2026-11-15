# Analisi della fattorizzazione LU: implementazione, ottimizzazioni e risultati

## 1. Obiettivo e sintesi

Il programma implementa una fattorizzazione LU senza pivoting e ne confronta
una versione seriale con diverse varianti OpenMP. Questa relazione descrive
il funzionamento del kernel, le ottimizzazioni incrementali e i risultati
presenti in [`lu_kernel_results.csv`](./lu_kernel_results.csv). Le misure
principali sono i tempi **del kernel**, non dell'intero processo.

Per `N=4000` la mediana dell'originale seriale è **60.087601 s**, quella della
versione con pragma OpenMP **38.538919 s**, quella SIMD **37.820548 s** e
quelle `static` e `static,1` rispettivamente **37.536679 s** e **37.481836 s**.
Il rapporto fra il tempo dell'originale e quello di `static` è **1.601×**.
Le differenze tra le ultime due versioni sono molto più piccole del guadagno
rispetto all'originale.

**Limite importante:** la matrice usata dal programma è singolare. Le misure
descrivono il comportamento di queste esecuzioni, ma non certificano la
correttezza numerica della fattorizzazione né le prestazioni su una matrice
generica non singolare.

## 2. Struttura del programma

Il programma opera su una matrice quadrata `A` di dimensione `N × N`,
normalmente di tipo `double`. Le dimensioni dei dataset sono definite in
[`lu.h`](./lu.h).

| Fase | Cosa avviene |
|---|---|
| Allocazione | Si riserva la matrice `A` tramite le macro PolyBench. |
| Inizializzazione | `init_array` imposta `A[i][j] = ((i+1)(j+1))/N`. |
| Misura | Le macro PolyBench avviano il timer immediatamente prima di `kernel_lu`. |
| Calcolo | `kernel_lu` modifica `A` sul posto, senza allocare matrici `L` e `U` separate. |
| Fine | Si ferma e stampa il timer; la stampa della matrice per evitare dead-code elimination dipende dalle macro PolyBench e dagli argomenti; infine si libera la memoria. |

### 2.1 Cosa fa `kernel_lu`

Il codice seriale è in [`lu.c`](./lu.c). L'idea è elaborare un pivot `k` alla
volta, normalizzare la parte destra della riga del pivot e aggiornare la
sottomatrice sottostante:

```text
per k = 0, ..., N-1:
    per j = k+1, ..., N-1:
        A[k][j] = A[k][j] / A[k][k]

    per i = k+1, ..., N-1:
        per j = k+1, ..., N-1:
            A[i][j] = A[i][j] - A[i][k] * A[k][j]
```

Il valore diagonale `A[k][k]` è il pivot. Le posizioni sopra la diagonale
contengono i coefficienti normalizzati di `U`; la diagonale di `U` è
implicitamente pari a uno. La parte inferiore, diagonale inclusa, contiene
`L`: questa è una convenzione di memorizzazione *in place* per `A = L U`.
Non si cercano pivot alternativi e non si scambiano righe.

Il ciclo su `k` è sequenziale: il pivot successivo dipende dagli aggiornamenti
precedenti. Per un dato `k`, invece, le colonne della riga pivot sono
indipendenti tra loro; una volta completata la normalizzazione, anche gli
aggiornamenti di righe diverse possono essere eseguiti in parallelo. Occorre
una barriera tra normalizzazione e aggiornamento e una seconda prima del pivot
successivo. Il lavoro aritmetico del kernel è **O(N³)** e la matrice occupa
**O(N²)** memoria; all'aumentare di `k` la sottomatrice si riduce.

## 3. Ottimizzazioni realizzate

Le versioni intermedie sono implementate in
[`lu_v2_steps.c`](./lu_v2_steps.c), selezionando `LU_STEP` da 1 a 5. Ogni
passaggio conserva quelli precedenti.

| Ordine | Nome nel CSV | Modifica | Obiettivo e limite |
|---:|---|---|---|
| 0 | `v0` | Originale seriale | Riferimento; il compilatore può comunque applicare proprie ottimizzazioni. |
| 1 | `v2` | Una regione `omp parallel` persistente e due `omp for` per pivot | Distribuire colonne e righe senza ricreare il team a ogni pivot. Le barriere implicite garantiscono le dipendenze, ma comportano attesa e sincronizzazione. |
| 2 | `v2_simd` | `omp for simd` per normalizzare e `omp simd` sul ciclo interno di aggiornamento | Esporre il parallelismo tra elementi alle istruzioni vettoriali. La direttiva non prova da sola che il compilatore abbia prodotto istruzioni SIMD. |
| 3 | `v2_recip` | Un reciproco `1/A[k][k]` per thread e pivot, poi moltiplicazioni | Ridurre il numero di divisioni; può modificare gli arrotondamenti e non garantisce un guadagno se il compilatore fa già la trasformazione. |
| 4 | `v2_hoist` | Caricamento di `A[i][k]` fuori dal ciclo interno su `j` | Evitare caricamenti ripetuti dello stesso valore per riga. |
| 5 | `v2_rows` | Puntatori `Ai` e `Ak` alle due righe | Rendere espliciti gli accessi contigui e potenzialmente semplificare l'indirizzamento. |
| 6 | `v2_restrict` | `restrict` sui puntatori alle righe | Comunicare al compilatore il vincolo di aliasing e favorire ottimizzazioni. Le righe di pivot e aggiornamento sono distinte perché `i > k`. |
| 7 | `v3`, **static** | Scheduling `static` esplicito, con diversa gestione degli indici | Dividere il lavoro in blocchi approssimativamente uguali; rispetto a `v2_restrict` non cambia esclusivamente la clausola di scheduling. |
| 8 | `v4`, **static,1** | `schedule(static,1)` sulle righe di aggiornamento | Assegnare le righe ciclicamente; normalizzazione ancora con `static`. |

`static` senza chunk esplicito distribuisce blocchi approssimativamente uguali;
`static,1` assegna una riga per volta in round-robin. Per esempio, con quattro
thread, quest'ultimo assegna al primo le righe 0, 4, 8, … e al secondo 1, 5,
9, …. L'assegnazione a blocchi può favorire la località, quella ciclica può
aiutare su lavoro irregolare. Qui, per un pivot fissato, le righe aggiornate
hanno cicli della stessa lunghezza: il beneficio di bilanciamento non è
garantito. I sorgenti delle due versioni sono [`lu_v3.c`](./lu_v3.c) e
[`lu_v4.c`](./lu_v4.c).

## 4. Metodologia di valutazione

| Aspetto | Misure disponibili o metodo |
|---|---|
| Dataset | `N=32`, 128, 1024, 2000, 4000. |
| Thread | `threads=4` in tutte le righe; `v0` rimane seriale. Il numero effettivo di thread del team non è registrato. |
| Ripetizioni | Tre per ogni versione e dimensione nei dati `perf`. Si riportano la mediana e, dove utile, minimo e massimo. |
| Copertura | 135 righe `perf` = 9 versioni × 5 dimensioni × 3 prove. Ci sono inoltre 81 righe Cachegrind, solo per `N=32`, 128 e 1024. Tutti gli `exit_status` riportati sono 0. |
| Tempo del kernel | `kernel_seconds` misura la chiamata `kernel_lu`, comprese creazione del team OpenMP e sincronizzazioni, ma non inizializzazione, stampa o deallocazione. |
| Tempo di processo | `wall_seconds` comprende l'intera invocazione e non è usato per gli speedup riportati qui. |
| Contatori | `perf` conta l'intero processo, compreso il runtime OpenMP; Cachegrind simula eventi dell'intero processo e altera l'esecuzione. I suoi tempi non sono tempi nativi. |

Lo script attuale [`perf_valgrind.sh`](./perf_valgrind.sh) abilita il timer
PolyBench e disabilita il flush cache preliminare di PolyBench. Il timer usa
`gettimeofday` e stampa sei decimali: non è un clock monotono e la risoluzione
limita soprattutto l'interpretazione dei dataset piccoli. Le misure del kernel
nelle righe `perf` sono comunque raccolte **sotto `perf`**, il cui overhead può
incidere. Il CSV non registra CPU, compilatore e flag effettivi, affinità,
frequenza, carico, data, ordine, numero effettivo di thread o multiplexing;
non è possibile ricostruire o verificare queste condizioni dal solo CSV.

Definiamo lo speedup per dimensione come:

```text
speedup(versione) = mediana(kernel_seconds originale) /
                    mediana(kernel_seconds versione)
```

Un valore superiore a 1 indica una mediana inferiore all'originale; è un
**rapporto di mediane**, non la mediana degli speedup per prova. I min–max
descrivono tre soli campioni, non intervalli di confidenza. Anche quando i
range non si sovrappongono, non costituiscono da soli una prova statistica.

## 5. Tempi e speedup

### 5.1 Mediane di tutte le versioni

Tempi in secondi, `perf`, 4 thread richiesti e tre ripetizioni. I valori sono
riportati con la precisione del timer:

| Versione | N=32 | N=128 | N=1024 | N=2000 | N=4000 |
|---|---:|---:|---:|---:|---:|
| Originale `v0` | 0.000031 | 0.001747 | 1.076367 | 7.419219 | 60.087601 |
| OpenMP `v2` | 0.000294 | 0.000939 | 0.597172 | 4.641734 | 38.538919 |
| SIMD | 0.000288 | 0.000884 | 0.582767 | 4.595963 | 37.820548 |
| + reciproco | 0.000292 | 0.000879 | 0.583056 | 4.574943 | 37.821681 |
| + hoisting | 0.000290 | 0.000838 | 0.568133 | 4.549707 | 37.554381 |
| + puntatori | 0.000291 | 0.000830 | 0.591544 | 4.545261 | 37.537228 |
| + `restrict` | 0.000287 | 0.000803 | 0.564595 | 4.544747 | 37.531454 |
| `static` | 0.000286 | 0.000819 | 0.575045 | 4.527864 | 37.536679 |
| `static,1` | 0.000299 | 0.001219 | 0.582838 | 4.522652 | 37.481836 |

Per `N=32`, l'originale impiega circa 31 microsecondi e le versioni parallele
circa 286–299 microsecondi: il lavoro è troppo piccolo per ammortizzare
creazione del team e sincronizzazioni. A `N=128` l'originale è più lento
delle versioni parallele mostrate, ma si tratta ancora di tempi inferiori a
pochi millisecondi: piccoli cambiamenti e precisione del timer pesano molto.

### 5.2 Confronto delle versioni principali

| N | Speedup OpenMP | Speedup SIMD | Speedup static | Speedup static,1 |
|---:|---:|---:|---:|---:|
| 32 | 0.105× | 0.108× | 0.108× | 0.104× |
| 128 | 1.860× | 1.976× | 2.133× | 1.433× |
| 1024 | 1.802× | 1.847× | 1.872× | 1.847× |
| 2000 | 1.598× | 1.614× | 1.639× | 1.640× |
| 4000 | 1.559× | 1.589× | 1.601× | 1.603× |

Ogni valore è calcolato come **tempo originale / tempo della versione**.

![Speedup per dimensione rispetto all'originale](./grafici_relazione_kernel/speedup_per_dataset.svg)

Per `N=4000`, il passaggio dall'originale a OpenMP riduce la mediana da
60.087601 a 38.538919 s. SIMD risulta più veloce della versione con soli
pragma nella campagna riportata: 37.820548 contro 38.538919 s, circa **1.86%**
in meno. Rispetto a SIMD, `static` ha una mediana inferiore di circa **0.75%**;
questa differenza aggrega più modifiche e non identifica il contributo di una
singola ottimizzazione.

### 5.3 Passaggi successivi e variabilità a N=4000

![Tempi mediani dei passaggi, N=4000](./grafici_relazione_kernel/tempi_passaggi_n4000.svg)

| Versione | Mediana kernel (s) | Min–max (s) | Speedup su originale |
|---|---:|---:|---:|
| Originale | 60.087601 | 60.007218–60.116156 | 1.000× |
| OpenMP | 38.538919 | 38.535083–38.546236 | 1.559× |
| SIMD | 37.820548 | 37.809507–37.846906 | 1.589× |
| + reciproco | 37.821681 | 37.806470–38.231522 | 1.589× |
| + hoisting | 37.554381 | 37.544441–37.559699 | 1.600× |
| + puntatori | 37.537228 | 37.532494–37.551365 | 1.601× |
| + `restrict` | 37.531454 | 37.524102–37.534687 | 1.601× |
| `static` | 37.536679 | 37.529436–37.556152 | 1.601× |
| `static,1` | 37.481836 | 37.480024–37.490305 | 1.603× |

Reciproco e SIMD hanno mediane quasi identiche; una delle prove del reciproco
è più lenta. Il successivo hoisting ha una mediana inferiore, mentre puntatori
e `restrict` apportano differenze ulteriormente piccole. Questi sono confronti
tra versioni cumulative: il compilatore può ottimizzare autonomamente, e
variabilità e ordine di esecuzione non consentono di attribuire con certezza
ogni differenza alla sola modifica indicata.

## 6. Confronto `static` e `static,1`

La differenza percentuale è `(T_static,1 / T_static − 1) × 100`: un valore
negativo indica una mediana più bassa per `static,1`.

| N | `static` mediana [min–max] s | `static,1` mediana [min–max] s | Differenza |
|---:|---:|---:|---:|
| 32 | 0.000286 [0.000286–0.000287] | 0.000299 [0.000296–0.000299] | +4.55% |
| 128 | 0.000819 [0.000811–0.000899] | 0.001219 [0.001208–0.001243] | +48.84% |
| 1024 | 0.575045 [0.553745–0.579151] | 0.582838 [0.561424–0.593323] | +1.36% |
| 2000 | 4.527864 [4.524003–4.533600] | 4.522652 [4.505185–4.523440] | −0.12% |
| 4000 | 37.536679 [37.529436–37.556152] | 37.481836 [37.480024–37.490305] | −0.15% |

![Differenza percentuale tra static,1 e static per dataset](./grafici_relazione_kernel/static_vs_static1_dataset.svg)

La differenza maggiore osservata è a `N=128`, ma riguarda tempi molto piccoli
e non dimostra che la stessa tendenza valga per altri input. Per `N=1024` i
range si sovrappongono; per `N=4000` le tre prove disponibili di `static,1`
sono inferiori a quelle di `static`, ma il divario tra mediane è solo
0.054843 s. Non si può stabilire una superiorità generale con tre prove e
senza controlli aggiuntivi su macchina, ordine, affinità e input.

## 7. Istruzioni, cicli e IPC

L'IPC (*instructions per cycle*) è calcolato per ogni prova come
`instructions / cycles`; la tabella e il grafico riportano la mediana dei tre
rapporti. Le colonne istruzioni e cicli sono invece mediane separate, espresse
in miliardi. Questi contatori riguardano il **processo completo**, non
soltanto il kernel, e aggregano anche l'attività dei thread e del runtime.

![IPC per le varie versioni, N=4000](./grafici_relazione_kernel/ipc_passaggi_n4000.svg)

| Versione, N=4000 | Istruzioni mediane (miliardi) | Cicli mediani (miliardi) | IPC mediano | Tempo kernel mediano (s) |
|---|---:|---:|---:|---:|
| Originale | 170.956 | 89.358 | 1.913 | 60.087601 |
| OpenMP | 174.394 | 226.936 | 0.769 | 38.538919 |
| SIMD | 110.533 | 222.516 | 0.497 | 37.820548 |
| + reciproco | 110.801 | 222.614 | 0.498 | 37.821681 |
| + hoisting | 89.650 | 220.970 | 0.406 | 37.554381 |
| + puntatori | 89.346 | 220.937 | 0.404 | 37.537228 |
| + `restrict` | 89.264 | 220.886 | 0.404 | 37.531454 |
| `static` | 89.379 | 220.843 | 0.405 | 37.536679 |
| `static,1` | 89.227 | 220.623 | 0.404 | 37.481836 |

L'originale ha IPC più alto, ma tempo trascorso maggiore: **IPC non è uno
speedup**. In parallelo i cicli sono aggregati tra thread e includono anche
runtime e sincronizzazione; più thread possono progredire contemporaneamente.
Inoltre, una singola istruzione vettoriale può eseguire più operazioni utili.

Tra SIMD e hoisting la mediana delle istruzioni passa da 110.533 a 89.650
miliardi, circa **−18.89%**, mentre i cicli passano da 222.516 a 220.970
miliardi. L'IPC scende da 0.497 a 0.406, anche se il tempo del kernel
diminuisce leggermente. I dati sono coerenti con meno istruzioni su cicli
simili, ma non isolano la causa né provano l'effettiva vettorizzazione.
Per `N=4000`, i cache miss mediani di processo sono circa 4.635 miliardi per
`static` e 4.700 per `static,1`: non bastano a spiegare una differenza di
tempo o a determinare da soli quale scheduling abbia migliore località.

## 8. Correttezza e considerazioni finali

La matrice inizializzata è il prodotto esterno dei vettori `i+1` e `j+1`,
scalato per `1/N`: per `N>1` ha rango uno ed è **singolare**. Senza pivoting,
il kernel può incontrare pivot nulli o numericamente instabili e produrre
`NaN` o `Inf`. Un `exit_status` pari a zero dice solo che il processo è
terminato senza errore segnalato: **non verifica che la fattorizzazione sia
corretta**. Anche differenze di arrotondamento tra divisione e reciproco
possono essere rilevanti.

In sintesi:

1. Nei dataset maggiori, le versioni OpenMP hanno tempi mediani inferiori
   all'originale; a `N=4000` lo speedup osservato di `static` è 1.601×.
2. Le micro-ottimizzazioni successive a SIMD mostrano miglioramenti di tempo
   piccoli e non uniformi, pur con una riduzione osservata delle istruzioni.
3. `static` e `static,1` sono molto vicini per `N=2000` e 4000; `static,1`
   ha mediane leggermente inferiori in queste misure. Non emerge una scelta
   universalmente migliore dai soli dati disponibili.
4. L'IPC, i contatori e la variazione su dataset piccoli vanno interpretati
   insieme al tempo del kernel e ai limiti della misura, non isolatamente.


