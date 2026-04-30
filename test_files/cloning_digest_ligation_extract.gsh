# Cloning macro fixture: Digest -> Ligation -> ExtractRegion
op {"Digest":{"input":"${seq_id}","enzymes":["BamHI"],"output_prefix":"${digest_prefix}"}}
op {"Ligation":{"inputs":["${digest_prefix}_1","${digest_prefix}_2"],"circularize_if_possible":false,"output_id":null,"protocol":"Sticky","output_prefix":"${ligation_prefix}","unique":false}}
op {"ExtractRegion":{"input":"${ligation_prefix}_1","from":${extract_from},"to":${extract_to},"output_id":"${output_id}"}}
