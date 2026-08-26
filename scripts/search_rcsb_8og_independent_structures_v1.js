const fs = require('fs');
const path = require('path');
const crypto = require('crypto');


const project = path.resolve(__dirname, '..');
const outDir = path.join(
  project,
  'data',
  'pipeline_workspace',
  '07_issue_resolution',
  'issue1_rcsb_search_2026-08-04'
);
const rawDir = path.join(outDir, 'raw_entry_metadata');
fs.mkdirSync(rawDir, { recursive: true });

const SEARCH_URL = 'https://search.rcsb.org/rcsbsearch/v2/query';
const DATA_URL = 'https://data.rcsb.org/rest/v1/core/entry/';


function sha256(value) {
  return crypto.createHash('sha256').update(value).digest('hex');
}


function csvEscape(value) {
  if (value === null || value === undefined) return '';
  const text = String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}


function terminal(attribute, operator, value) {
  return { type: 'terminal', service: 'text', parameters: { attribute, operator, value } };
}


function makeQuery(nodes) {
  return {
    query: nodes.length === 1 ? nodes[0] : { type: 'group', logical_operator: 'and', nodes },
    request_options: {
      paginate: { start: 0, rows: 10000 },
      results_content_type: ['experimental'],
      results_verbosity: 'minimal'
    },
    return_type: 'entry'
  };
}


const eightOg = terminal(
  'rcsb_polymer_entity_container_identifiers.chem_comp_monomers',
  'exact_match',
  '8OG'
);
const proteinCountZero = terminal('rcsb_entry_info.polymer_entity_count_protein', 'equals', 0);
const xray = terminal('exptl.method', 'exact_match', 'X-RAY DIFFRACTION');

const queries = {
  all_polymer_8OG: makeQuery([eightOg]),
  protein_free_polymer_8OG: makeQuery([eightOg, proteinCountZero]),
  protein_free_xray_polymer_8OG: makeQuery([eightOg, proteinCountZero, xray])
};


async function fetchText(url, options = {}) {
  const started = new Date().toISOString();
  const response = await fetch(url, options);
  const text = await response.text();
  if (!response.ok) {
    throw new Error(`${url}: HTTP ${response.status}: ${text.slice(0, 500)}`);
  }
  return {
    text,
    status: response.status,
    contentType: response.headers.get('content-type') || '',
    etag: response.headers.get('etag') || '',
    lastModified: response.headers.get('last-modified') || '',
    retrievedAtUtc: started
  };
}


function getIdentifiers(parsed) {
  return (parsed.result_set || []).map(row => row.identifier).sort();
}


async function main() {
  const searchManifest = [];
  const searchResults = {};

  for (const [name, query] of Object.entries(queries)) {
    const requestText = JSON.stringify(query, null, 2) + '\n';
    const requestPath = path.join(outDir, `${name}_request.json`);
    fs.writeFileSync(requestPath, requestText, 'utf8');

    const response = await fetchText(SEARCH_URL, {
      method: 'POST',
      headers: { 'content-type': 'application/json' },
      body: JSON.stringify(query)
    });
    const parsed = JSON.parse(response.text);
    const normalizedResponse = JSON.stringify(parsed, null, 2) + '\n';
    const responsePath = path.join(outDir, `${name}_response.json`);
    fs.writeFileSync(responsePath, normalizedResponse, 'utf8');
    const identifiers = getIdentifiers(parsed);
    searchResults[name] = identifiers;
    searchManifest.push({
      query_name: name,
      endpoint: SEARCH_URL,
      retrieved_at_utc: response.retrievedAtUtc,
      http_status: response.status,
      content_type: response.contentType,
      etag: response.etag,
      last_modified: response.lastModified,
      total_count: parsed.total_count,
      response_sha256: sha256(Buffer.from(normalizedResponse, 'utf8')),
      identifiers_json: JSON.stringify(identifiers)
    });
  }

  const proteinFreeIds = searchResults.protein_free_polymer_8OG;
  const metadataRows = [];
  for (const pdbId of proteinFreeIds) {
    const response = await fetchText(`${DATA_URL}${pdbId}`);
    const parsed = JSON.parse(response.text);
    const normalized = JSON.stringify(parsed, null, 2) + '\n';
    const rawPath = path.join(rawDir, `${pdbId}.json`);
    fs.writeFileSync(rawPath, normalized, 'utf8');

    const info = parsed.rcsb_entry_info || {};
    const accession = parsed.rcsb_accession_info || {};
    const citations = parsed.citation || [];
    const primaryCitation = citations.find(row => row.id === 'primary') || citations[0] || {};
    const methods = (parsed.exptl || []).map(row => row.method).filter(Boolean);
    const resolutions = info.resolution_combined || [];
    metadataRows.push({
      pdb_id: pdbId,
      title: parsed.struct?.title || '',
      methods: methods.join(' | '),
      resolution_A: resolutions.join(' | '),
      initial_release_date: accession.initial_release_date || '',
      protein_entity_count: info.polymer_entity_count_protein ?? '',
      DNA_entity_count: info.polymer_entity_count_DNA ?? '',
      RNA_entity_count: info.polymer_entity_count_RNA ?? '',
      nucleic_acid_entity_count: info.polymer_entity_count_nucleic_acid ?? '',
      polymer_entity_count: info.polymer_entity_count ?? '',
      primary_citation_doi: primaryCitation.pdbx_database_id_DOI || '',
      rcsb_entry_url: `https://www.rcsb.org/structure/${pdbId}`,
      data_api_url: `${DATA_URL}${pdbId}`,
      metadata_retrieved_at_utc: response.retrievedAtUtc,
      metadata_sha256: sha256(Buffer.from(normalized, 'utf8'))
    });
  }

  const manifestHeaders = Object.keys(searchManifest[0]);
  fs.writeFileSync(
    path.join(outDir, 'search_manifest.csv'),
    '\ufeff' + [
      manifestHeaders.join(','),
      ...searchManifest.map(row => manifestHeaders.map(header => csvEscape(row[header])).join(','))
    ].join('\n') + '\n',
    'utf8'
  );

  const metadataHeaders = Object.keys(metadataRows[0]);
  fs.writeFileSync(
    path.join(outDir, 'protein_free_8OG_entry_metadata.csv'),
    '\ufeff' + [
      metadataHeaders.join(','),
      ...metadataRows.map(row => metadataHeaders.map(header => csvEscape(row[header])).join(','))
    ].join('\n') + '\n',
    'utf8'
  );

  const xrayIds = new Set(searchResults.protein_free_xray_polymer_8OG);
  const xrayRows = metadataRows.filter(row => xrayIds.has(row.pdb_id));
  const audit = {
    status: 'PASS_RCSB_SEARCH_SNAPSHOT',
    as_of_utc: new Date().toISOString(),
    official_search_endpoint: SEARCH_URL,
    query_counts: Object.fromEntries(searchManifest.map(row => [row.query_name, row.total_count])),
    protein_free_xray_ids: [...xrayIds].sort(),
    protein_free_xray_titles: Object.fromEntries(xrayRows.map(row => [row.pdb_id, row.title])),
    interpretation_gate: {
      automated_search_complete_for_defined_query: true,
      title_and_method_metadata_saved: true,
      manual_duplex_and_pair_partner_review_required: true,
      protein_bound_entries_are_not_independent_free_DNA_replicates: true,
      NMR_entries_require_separate_stratum: true
    }
  };
  fs.writeFileSync(
    path.join(outDir, 'rcsb_8OG_search_summary.json'),
    JSON.stringify(audit, null, 2) + '\n',
    'utf8'
  );
  process.stdout.write(JSON.stringify(audit, null, 2) + '\n');
}


main().catch(error => {
  console.error(error);
  process.exit(1);
});
