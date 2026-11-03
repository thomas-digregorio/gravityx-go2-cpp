import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const workbookPath = "C:/Users/thoma/Documents/goc2-ac-score-check/final.xlsx";
const input = await FileBlob.load(workbookPath);
const workbook = await SpreadsheetFile.importXlsx(input);

const sheets = await workbook.inspect({
  kind: "sheet",
  include: "id,name",
  maxChars: 12000,
});
console.log("SHEETS");
console.log(sheets.ndjson);

const hits = await workbook.inspect({
  kind: "match",
  searchTerm: "C2FEN19402|19402|Division 1|Div 1|MSpp|459401|549955|953207|395220|444334",
  options: { useRegex: true, maxResults: 500 },
  maxChars: 30000,
  summary: "GO2 19,402-bus score references",
});
console.log("MATCHES");
console.log(hits.ndjson);
