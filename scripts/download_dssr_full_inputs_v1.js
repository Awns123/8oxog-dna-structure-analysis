const crypto = require("crypto");
const fs = require("fs");
const path = require("path");

const outputDir = path.resolve(
  __dirname,
  "..",
  "data",
  "pipeline_workspace",
  "01_raw_mmcif",
  "full_v1_2026-08-04",
);

const referenceIds = [
  "119D", "158D", "167D", "1BNA", "1D29", "1D49", "1D65",
  "1D89", "1D98", "1DN9", "1JGR", "2BNA", "3BSE", "3IXN",
  "463D", "476D", "477D", "4C64",
];
const targetIds = ["111D", "178D", "183D"];
const inputs = [
  ...referenceIds.map((id) => ({
    pdbId: id,
    role: "normal_reference_v1",
    name: `${id}.cif`,
    url: `https://files.rcsb.org/download/${id}.cif`,
  })),
  ...targetIds.map((id) => ({
    pdbId: id,
    role: id === "111D" ? "unoxidized_GA_matched_analog"
      : id === "178D" ? "8OG_A_case" : "8OG_C_ASU_audit",
    name: `${id}.cif`,
    url: `https://files.rcsb.org/download/${id}.cif`,
  })),
  {
    pdbId: "183D",
    role: "8OG_C_primary_assembly1",
    name: "183D-assembly1.cif",
    url: "https://files.rcsb.org/download/183D-assembly1.cif",
  },
];

function csvCell(value) {
  const text = value == null ? "" : String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function sha256(bytes) {
  return crypto.createHash("sha256").update(bytes).digest("hex");
}

function verifyLockedManifest(manifestPath, rows) {
  const lines = fs.readFileSync(manifestPath, "utf8").trim().split(/\r?\n/);
  const headers = lines[0].split(",");
  const filenameIndex = headers.indexOf("filename");
  const sizeIndex = headers.indexOf("size_bytes");
  const shaIndex = headers.indexOf("sha256");
  if ([filenameIndex, sizeIndex, shaIndex].some((index) => index < 0)) {
    throw new Error("locked input manifest is missing filename, size_bytes, or sha256");
  }
  const locked = new Map(
    lines.slice(1).map((line) => {
      const cells = line.split(",");
      return [
        cells[filenameIndex],
        { size: Number(cells[sizeIndex]), sha: cells[shaIndex].toLowerCase() },
      ];
    }),
  );
  for (const row of rows) {
    const expected = locked.get(row.filename);
    if (!expected) {
      throw new Error(row.filename + ": missing from locked input manifest");
    }
    if (expected.size !== row.size_bytes || expected.sha !== row.sha256.toLowerCase()) {
      throw new Error(
        row.filename + ": downloaded bytes do not match locked manifest " +
        "(size " + row.size_bytes + "/" + expected.size +
        ", sha " + row.sha256 + "/" + expected.sha + ")",
      );
    }
  }
  if (locked.size !== rows.length) {
    throw new Error(
      "locked manifest row count " + locked.size +
      " != downloaded row count " + rows.length,
    );
  }
}

async function downloadOne(item) {
  const finalPath = path.join(outputDir, item.name);
  const partialPath = `${finalPath}.partial`;
  if (fs.existsSync(partialPath)) {
    throw new Error(item.name + ": partial file remains; inspect or remove it before retrying");
  }
  if (fs.existsSync(finalPath)) {
    const bytes = fs.readFileSync(finalPath);
    const prefix = bytes.subarray(0, 32).toString("utf8");
    if (!prefix.startsWith("data_")) {
      throw new Error(item.name + ": existing file is not an uncompressed mmCIF");
    }
    return {
      pdb_id: item.pdbId,
      filename: item.name,
      role: item.role,
      source_url: item.url,
      retrieved_at_utc: fs.statSync(finalPath).mtime.toISOString(),
      http_status: "REUSED_EXISTING",
      content_type: "chemical/x-cif",
      content_length_header: bytes.length,
      etag: "",
      last_modified: "",
      size_bytes: bytes.length,
      sha256: sha256(bytes),
    };
  }

  const startedAt = new Date().toISOString();
  const response = await fetch(item.url, {
    headers: { "User-Agent": "8oxog-dna-structure-analysis/1.0" },
  });
  if (!response.ok) {
    throw new Error(`${item.name}: HTTP ${response.status}`);
  }
  const bytes = Buffer.from(await response.arrayBuffer());
  const prefix = bytes.subarray(0, 32).toString("utf8");
  if (!prefix.startsWith("data_")) {
    throw new Error(`${item.name}: response is not an uncompressed mmCIF (prefix=${JSON.stringify(prefix)})`);
  }
  const declaredLength = response.headers.get("content-length");
  if (declaredLength && Number(declaredLength) !== bytes.length) {
    throw new Error(`${item.name}: content-length mismatch ${declaredLength} != ${bytes.length}`);
  }

  fs.writeFileSync(partialPath, bytes, { flag: "wx" });
  fs.renameSync(partialPath, finalPath);

  return {
    pdb_id: item.pdbId,
    filename: item.name,
    role: item.role,
    source_url: item.url,
    retrieved_at_utc: startedAt,
    http_status: response.status,
    content_type: response.headers.get("content-type") || "",
    content_length_header: declaredLength || "",
    etag: response.headers.get("etag") || "",
    last_modified: response.headers.get("last-modified") || "",
    size_bytes: bytes.length,
    sha256: sha256(bytes),
  };
}

async function main() {
  fs.mkdirSync(outputDir, { recursive: true });
  const rows = [];
  for (const item of inputs) {
    const row = await downloadOne(item);
    rows.push(row);
    process.stdout.write(`${row.filename}\t${row.size_bytes}\t${row.sha256}\n`);
  }

  const headers = Object.keys(rows[0]);
  const csv = [
    headers.join(","),
    ...rows.map((row) => headers.map((key) => csvCell(row[key])).join(",")),
  ].join("\n") + "\n";
  const manifestPath = path.join(outputDir, "input_manifest_full_v1.csv");
  if (fs.existsSync(manifestPath)) {
    verifyLockedManifest(manifestPath, rows);
    const logDir = path.resolve(outputDir, "..", "..", "08_logs", "downloads");
    fs.mkdirSync(logDir, { recursive: true });
    fs.writeFileSync(path.join(logDir, "input_manifest_rerun_full_v1.csv"), csv, "utf8");
    process.stdout.write("LOCKED_MANIFEST_MATCH\t22/22\n");
  } else {
    fs.writeFileSync(manifestPath, csv, { encoding: "utf8", flag: "wx" });
  }
}

main().catch((error) => {
  process.stderr.write(`${error.stack || error}\n`);
  process.exitCode = 1;
});
