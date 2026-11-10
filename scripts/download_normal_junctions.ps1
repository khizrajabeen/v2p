# Download normal-tissue splice junction data from GTEx and recount3.
# Windows PowerShell 5.1. Both steps resume if interrupted.

$Root = "D:\v2p\incoming\normal_junctions"
New-Item -ItemType Directory -Force "$Root\gtex"     | Out-Null
New-Item -ItemType Directory -Force "$Root\recount3" | Out-Null

# ---------------------------------------------------------------- GTEx v8
# One matrix: junction read counts, all samples, all tissues. 4.25 GB.
$GtexUrl  = "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/" +
            "GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
$GtexFile = "$Root\gtex\GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"

Write-Host "GTEx v8 junctions -> $GtexFile"
curl.exe -fL --retry 3 -C - -o $GtexFile $GtexUrl

# Sample attributes, to map each sample column to its tissue. 11 MB.
curl.exe -fL --retry 3 -C - -o "$Root\gtex\GTEx_v8_SampleAttributesDS.txt" `
  "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

# ------------------------------------------------------------- recount3
# 32 tissues. Each has an .RR.gz (junction coordinates, ~32 MB) and an
# .MM.gz (per-sample counts, 0.3-2.8 GB). Only RR is fetched here: the
# full MM set is 37 GB and is only needed once a support threshold is set.
$Bucket = "https://recount-opendata.s3.amazonaws.com/"
$Prefix = "recount3/release/human/data_sources/gtex/junctions/"

$keys  = @()
$token = $null
do {
    $url = "$Bucket" + "?list-type=2&prefix=$Prefix&max-keys=1000"
    if ($token) { $url += "&continuation-token=" + [uri]::EscapeDataString($token) }
    [xml]$xml = (Invoke-WebRequest -Uri $url -UseBasicParsing).Content
    $keys += $xml.ListBucketResult.Contents | Where-Object { $_.Key -like "*.RR.gz" }
    if ($xml.ListBucketResult.IsTruncated -eq "true") {
        $token = $xml.ListBucketResult.NextContinuationToken
    } else {
        $token = $null
    }
} while ($token)

$totalGb = ($keys | Measure-Object -Property Size -Sum).Sum / 1GB
Write-Host ("recount3: {0} junction files, {1:N2} GB" -f $keys.Count, $totalGb)

foreach ($k in $keys) {
    $name = Split-Path $k.Key -Leaf
    $dest = Join-Path "$Root\recount3" $name
    if ((Test-Path $dest) -and ((Get-Item $dest).Length -eq [int64]$k.Size)) {
        Write-Host "have $name"; continue
    }
    Write-Host ("{0}  {1:N0} MB" -f $name, ([int64]$k.Size / 1MB))
    curl.exe -fL --retry 3 -C - -o $dest ($Bucket + $k.Key)
}

Write-Host "Done. Files are under $Root"
