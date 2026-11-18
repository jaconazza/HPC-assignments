// Generate the two comparison charts with Node.js, without external packages.
const fs = require("node:fs");
const path = require("node:path");

const datasets = [
  { name: "MINI", csv: "MINI_DATASET", n: 32 },
  { name: "SMALL", csv: "SMALL_DATASET", n: 128 },
  { name: "STANDARD", csv: "STANDARD_DATASET", n: 1024 },
  { name: "LARGE", csv: "LARGE_DATASET", n: 2000 },
  { name: "EXTRALARGE", csv: "EXTRALARGE_DATASET", n: 4000 },
];
const blocks = [8, 16, 32, 64, 128];

function readCsv(filename) {
  const lines = fs.readFileSync(path.join(__dirname, "..", filename), "utf8")
    .trim().split(/\r?\n/);
  const headers = lines.shift().split(",");
  return lines.map((line) => {
    const values = line.split(",");
    if (values.length !== headers.length) {
      throw new Error(`Riga CSV non valida in ${filename}: ${line}`);
    }
    return Object.fromEntries(headers.map((header, i) => [header, values[i]]));
  });
}

function median(rows, field, description) {
  if (!rows.length) throw new Error(`Misure mancanti: ${description}`);
  const values = rows.map((row) => {
    const value = Number(row[field]);
    if (row.exit_status !== "0" || row[field] === "" || !Number.isFinite(value) || value <= 0) {
      throw new Error(`Misura non valida: ${description}`);
    }
    return value;
  }).sort((a, b) => a - b);
  const middle = Math.floor(values.length / 2);
  return values.length % 2 ? values[middle] : (values[middle - 1] + values[middle]) / 2;
}

const kernel = readCsv("lu_kernel_results.csv");
const tiled = readCsv("lu_tile_results.csv");
const results = datasets.map((dataset) => {
  const kernelMedian = (version) => median(
    kernel.filter((row) => row.tool === "perf" && row.threads === "4" &&
      Number(row.n) === dataset.n && row.version === version),
    "kernel_seconds", `${dataset.name}, ${version}`
  );
  const candidates = blocks.map((block) => ({
    block,
    time: median(tiled.filter((row) => row.tool === "perf" && row.nthreads === "4" &&
      row.dataset === dataset.csv && Number(row.block) === block),
    "exec_time_s", `${dataset.name}, blocco ${block}`),
  }));
  candidates.sort((a, b) => a.time - b.time || a.block - b.block);
  return {
    ...dataset,
    baseline: kernelMedian("v0"),
    simd: kernelMedian("v2_simd"),
    static: kernelMedian("v3"),
    tiled: candidates[0].time,
    block: candidates[0].block,
  };
});

const series = [
  { key: "simd", label: "v2_simd", color: "#16AFCB" },
  { key: "static", label: "v3 - schedule(static)", color: "#E94F85" },
  { key: "tiled", label: "Versione a blocchi", color: "#48A967" },
];

function escapeXml(text) {
  return String(text).replace(/&/g, "&amp;").replace(/</g, "&lt;")
    .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}

function chart({ title, yLabel, ticks, position, value, format, filename, reference, chartSeries = series }) {
  const width = 1100, height = 700;
  const left = 115, right = 985, top = 145, bottom = 470;
  const x = (index) => left + index * (right - left) / (results.length - 1);
  const y = (v) => bottom - position(v) * (bottom - top);
  const parts = [
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="title desc">`,
    `<title id="title">${escapeXml(title)}</title>`,
    `<desc id="desc">Confronto delle mediane dei tempi del kernel per cinque dataset e ${chartSeries.length} versioni, con 4 thread richiesti. La versione a blocchi usa il blocco con la mediana più bassa per ciascun dataset.</desc>`,
    `<defs>`,
    `<linearGradient id="brand" x1="0%" y1="0%" x2="100%" y2="0%">`,
    `<stop offset="0%" stop-color="#EF4C91"/>`,
    `<stop offset="55%" stop-color="#8B4DCE"/>`,
    `<stop offset="100%" stop-color="#3C56E8"/>`,
    `</linearGradient>`,
    `<linearGradient id="titlePanel" x1="0%" y1="0%" x2="100%" y2="100%">`,
    `<stop offset="0%" stop-color="#111936"/>`,
    `<stop offset="100%" stop-color="#27245F"/>`,
    `</linearGradient>`,
    `</defs>`,
    `<rect width="100%" height="100%" fill="#fff"/>`,
    `<rect x="4" y="4" width="1092" height="692" rx="18" fill="none" stroke="url(#brand)" stroke-width="7"/>`,
    `<g font-family="Arial, Helvetica, sans-serif" fill="#192344">`,
    `<rect x="28" y="22" width="224" height="30" rx="15" fill="url(#brand)"/>`,
    `<text x="140" y="42" text-anchor="middle" font-size="13" font-weight="bold" letter-spacing="1.1" fill="#fff">LU DECOMPOSITION</text>`,
    `<rect x="247" y="20" width="606" height="48" rx="24" fill="url(#titlePanel)"/>`,
    `<text x="550" y="51" text-anchor="middle" font-size="23" font-weight="bold" fill="#fff">${escapeXml(title)}</text>`,
    `<text x="${width / 2}" y="99" text-anchor="middle" font-size="14" fill="#56617D">Misure perf · 4 thread richiesti · mediana delle ripetizioni</text>`,
  ];
  for (const tick of ticks) {
    parts.push(`<line x1="${left}" y1="${y(tick)}" x2="${right}" y2="${y(tick)}" stroke="#DDE2EC" stroke-width="1"/>`);
    parts.push(`<text x="${left - 12}" y="${y(tick) + 5}" text-anchor="end" font-size="13" fill="#49546F">${escapeXml(format(tick))}</text>`);
  }
  if (reference !== undefined) {
    parts.push(`<line x1="${left}" y1="${y(reference)}" x2="${right}" y2="${y(reference)}" stroke="#A34CCB" stroke-width="2" stroke-dasharray="7 6"/>`);
  }
  parts.push(`<path d="M ${left} ${top} V ${bottom} H ${right}" fill="none" stroke="#35405F" stroke-width="1.5"/>`);
  parts.push(`<text x="28" y="${(top + bottom) / 2}" text-anchor="middle" font-size="14" fill="#35405F" transform="rotate(-90 28 ${(top + bottom) / 2})">${escapeXml(yLabel)}</text>`);
  results.forEach((result, i) => {
    parts.push(`<text x="${x(i)}" y="${bottom + 24}" text-anchor="middle" font-size="13" font-weight="bold" fill="#202B4C">${result.name}</text>`);
    parts.push(`<text x="${x(i)}" y="${bottom + 43}" text-anchor="middle" font-size="12" fill="#49546F">N=${result.n}</text>`);
    parts.push(`<text x="${x(i)}" y="${bottom + 61}" text-anchor="middle" font-size="12" fill="#68728A">blocco ${result.block}</text>`);
  });
  for (const item of chartSeries) {
    const points = results.map((result, i) => `${x(i)},${y(value(result, item.key))}`).join(" ");
    parts.push(`<polyline points="${points}" fill="none" stroke="${item.color}" stroke-width="3.5" stroke-linecap="round" stroke-linejoin="round"/>`);
    results.forEach((result, i) => {
      const v = value(result, item.key);
      parts.push(`<circle cx="${x(i)}" cy="${y(v)}" r="6" fill="${item.color}" stroke="#fff" stroke-width="2"><title>${escapeXml(`${item.label}, ${result.name}: ${v.toPrecision(6)}${item.key === "tiled" ? ` (blocco ${result.block})` : ""}`)}</title></circle>`);
    });
  }
  const legendWidth = chartSeries.length > 3 ? 900 : 730;
  const legendLeft = (width - legendWidth) / 2;
  const legendStep = legendWidth / chartSeries.length;
  parts.push(`<rect x="${legendLeft}" y="568" width="${legendWidth}" height="46" rx="23" fill="#F6F7FB" stroke="#E2E5EF"/>`);
  chartSeries.forEach((item, i) => {
    const legendX = chartSeries.length > 3 ? legendLeft + 15 + i * legendStep : 215 + i * 245;
    parts.push(`<line x1="${legendX}" y1="591" x2="${legendX + 34}" y2="591" stroke="${item.color}" stroke-width="3.5" stroke-linecap="round"/>`);
    parts.push(`<circle cx="${legendX + 17}" cy="591" r="5" fill="${item.color}" stroke="#fff" stroke-width="1.5"/>`);
    parts.push(`<text x="${legendX + 44}" y="596" font-size="13" font-weight="bold" fill="#303B5B">${escapeXml(item.label)}</text>`);
  });
  parts.push(`<text x="${width / 2}" y="642" text-anchor="middle" font-size="12" fill="#68728A">Blocchi scelti sugli stessi dati mostrati; misure provenienti da due raccolte distinte.</text>`);
  parts.push(`<text x="${width / 2}" y="662" text-anchor="middle" font-size="12" fill="#68728A">${escapeXml(reference === undefined ? "Tempo del kernel in secondi, scala logaritmica." : "Speed-up = mediana v0 / mediana versione; la linea tratteggiata indica 1×.")}</text>`);
  parts.push("</g></svg>");
  fs.writeFileSync(path.join(__dirname, filename), parts.join("\n") + "\n");
}

chart({
  title: "Confronto dei tempi del kernel LU",
  yLabel: "Tempo mediano del kernel (s, scala logaritmica)",
  ticks: [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100],
  position: (v) => (Math.log10(v) + 5) / 7,
  value: (result, key) => result[key],
  format: (v) => v.toString(),
  filename: "confronto_tempi_kernel_lu.svg",
  chartSeries: [
    { key: "baseline", label: "v0 (sequenziale)", color: "#808080" },
    ...series,
  ],
});
chart({
  title: "Confronto degli speed-up LU",
  yLabel: "Speed-up rispetto a v0",
  ticks: [0, 1, 2, 3, 4, 5],
  position: (v) => v / 5,
  value: (result, key) => result.baseline / result[key],
  format: (v) => `${v}×`,
  filename: "confronto_speedup_lu.svg",
  reference: 1,
});
