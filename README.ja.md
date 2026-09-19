# UniScFlow: one-command, reliable automated reconstruction and remapping of public single-cell and single-nucleus RNA-seq datasets from raw reads

<p align="center">
  <img src="docs/assets/uniscflow-emblem.svg" alt="UniScFlow emblem: platform-aware public scRNA-seq reprocessing" width="96%">
</p>

[English README](README.md)

**公開accessionからcount matrixへ。raw readの解釈もワークフローに組み込みます。**

`UniScFlow` は、GEO/SRA/ENAに公開されたsingle-cell and single-nucleus RNA-seq（sc-/snRNA-seq）をraw readから再解析するツールです。repository metadataと、raw readの構造・配列を直接調べた証拠を組み合わせ、platform、barcode/UMI geometry、logical read roleを推定します。PRJNAまたはGSE accessionと適切なreferenceを指定すると、`uniscflow --mode all` がinput取得、mapper-ready inputの再構成、対応経路のmappingまでをワンコマンドで実行します。

強みは、downloadやaligner実行に加えて、**公開されたread自体の役割を解釈すること**です。数字付きFASTQ suffixがbarcode、UMI、cDNA、indexの役割に揃っていなくても、利用可能な証拠を照合し、割り当てを記録してSTARsoloまたはSTAR + featureCounts用inputを準備します。chemistry定義とbarcode whitelistを渡すと、project内のsample別layoutを含む10x chemistry-aware検証も行えます。

**10x Genomics、Drop-seq、Seq-Well、Smart-seq2** はautomatic mappingに対応しています。[対応platform表](#対応platform)には、追加resourceを必要とする21 profileと、ユーザーがgeometryを明示するgeneric droplet-UMI routeも示しています。自動mappingに必要な証拠やresourceが不足する場合は、blockerと次の手順を報告します。Cell Ranger mappingはoptional legacy supportであり、default workflowではありません。

## まずここから

| 目的 | 推奨する進め方 |
| --- | --- |
| installを確認する | environmentを作成してUniScFlowをinstallし、`uniscflow --help` と `uniscflow --version` を実行します。full run前には、予定しているものと同じpath・設定で `--mode check` を使います。 |
| 軽量public demoを試す | [Smart-seq2版（英語）](docs/tutorials/quickstart_smartseq2_prjna701252.md)（3 GSM、sample map不要）か[10x版（英語）](docs/tutorials/quickstart_10x_prjna825585.md)（1 GSM、chemistry定義とbarcode whitelistを使用）を選びます。どちらにもDocker版があります。 |
| 自分のprojectを実行する | 対応するSTAR indexとGTFを用意し、SRA・展開FASTQ・mapper outputを同時に置ける容量を確保し、まず `--sample-alias` で確認済みsubsetを実行します。public 10x dataではCell Ranger chemistry定義とbarcode fileを強く推奨します。 |
| 終了結果を確認する | `.uniscflow_halt_after_download.json`、`mapper_inputs_manifest.tsv`、`mapper_run_manifest.tsv` を確認します。exit 0は、検証済みmapped outputだけでなく、意図したdocumented halt (recognized stop)やunsupported stopを表すこともあります。 |

以下のcommand例にある `/path/to/...` はplaceholderです。すべて実在pathに置き換え、段階実行を再開するときはmetadata、download script、raw data、mapperの各directoryを前段と同じにしてください。

## ワンコマンドworkflow

対応済みplatformでは、`--mode all` が基本のworkflowです。

<p align="center">
  <img src="docs/assets/uniscflow-one-command-workflow.svg" alt="UniScFlow one-command workflow: public IDからplatform推定、STAR-based mapping、安全なhalt、audit-ready outputまで" width="96%">
</p>

1つのcommandでsample/run関係を解決し、raw inputを取得・検証し、その構造を推定して対象sampleを指定referenceへmappingします。manifestには、選択accession、input file、推定read role、mapper command、outputの対応を残します。外部manifestやvendor固有preprocessingが必要な場合は、不足情報と次の手順を報告します。

下の詳細workflow mapは、`--mode all` が内部でmetadata取得、raw input downloadとintegrity check、platform/read構造推定、canonical mapper input生成、STARベースmappingまたは理由付きhalt、監査可能なreport作成へ展開される流れを示しています。

<p align="center">
  <img src="docs/assets/uniscflow-workflow.svg" alt="UniScFlowの詳細workflow map: metadata取得、raw input integrity check、platform/read構造推定、canonical mapper input生成、STAR-based mappingまたは理由付きhalt、output report作成" width="96%">
</p>

すべてのinputが全stageを通るわけではありません。FASTQ経路はcanonical linkを使い、raw barcode/UMI tagを検証済みの適格なBAMはFASTQ変換せずSAM streamとしてSTARsoloへ渡します。recognized stopはmapper preparation前にも起こり、matrixではなくreportを出力します。その他のmodeは確認や再開のために個々のstageを実行する補助modeです。通常は `--mode all` から始めます。

## Tutorial

まず **Lightweight public demo** から始めてください。Smart-seq2版と10x版は選択肢であり、両方を順番に実行する必要はありません。

| Version | 公開inputとmapping | Docker |
| --- | --- | --- |
| [Smart-seq2版（英語）](docs/tutorials/quickstart_smartseq2_prjna701252.md) | PRJNA701252の3 GSM（各1細胞）。platformとcell granularityを自動判定し、STAR + featureCountsで行列を作成します。sample map・barcode whitelistは不要です。 | [Smart-seq2 Docker版（英語）](docs/tutorials/docker_smartseq2_prjna701252.md) |
| [10x版（英語）](docs/tutorials/quickstart_10x_prjna825585.md) | PRJNA825585の1 GSMに属する全12 run。platform・chemistry・read roleを自動判定し、STARsoloで行列を作成します。chemistry定義とwhitelistの準備方法も説明しています。 | [10x Docker版（英語）](docs/tutorials/docker_10x_prjna825585.md) |

どちらもmouse STAR indexと対応GTFが必要です。「軽量」は選択する公開inputの規模を指し、whole-genome STAR indexのメモリ要件が小さくなるわけではありません。[tutorial guide（英語）](docs/tutorials/README.md)で共通の前提条件と出力構成を確認できます。

その他のworkflow:

| Tutorial | 内容 |
| --- | --- |
| [bulk + single-cell混在demo（英語）](docs/tutorials/mixed_bulk_gex_prjna949947.md) | PRJNA949947のbulk 2 GSMとSmart-seq2の1細胞をまとめて指定（約0.50 GB）。bulkを識別して対象外とし、single-cell GEXだけをmappingします。[Docker版（英語）](docs/tutorials/docker_mixed_bulk_gex_prjna949947.md)。 |
| [documented halt (recognized stop) demo（英語）](docs/tutorials/documented_halt.md) | vendor manifestや実験固有resourceがないと安全にmappingできないplatformで、download後に理由付きで停止する挙動を確認します。 |
| [複数projectのmouse tutorial（英語）](docs/examples/mouse_validation_tutorial.md) | PRJNA701252のmixed-platform例を含む、複数のpublic mouse PRJNAを処理します。実行scriptは `docs/examples/mouse_validation_tutorial.sh` です。 |
| [Tabula Muris Senis Brain Non-Myeloid tutorial（英語）](docs/examples/tabula_muris_senis_brain_nonmyeloid_tutorial.md) | GEO/SRA外のpublic FASTQをS3 directory構造を保って取得し、mapping layoutを作る例です。 |

`docs/tutorials/` はinstall確認と小規模な初回実行、`docs/examples/` はより大きな処理を説明します。まずいずれかの軽量demoを実行し、大きな例に進む前に必要なstorageとreferenceを確認してください。

## このworkflowが解決したい問題

複数の公共sc-/snRNA-seq datasetを均一に解析・比較・統合したい場合、GEOなどに上がっているcount matrixをそのまま使うだけでは不十分なことがあります。reference genome、gene annotation、gene ID、cell calling、filtering、mapperの違いが混ざってしまうためです。それぞれのmatrixは正しくても、dataset間で遺伝子定義や処理ruleが揃っていなければ、統合解析では大きな問題になります。

raw readを共通referenceへremappingすると、この処理を標準化できます。その際、alignmentの前に手間がかかるのが、sample/run関係の復元、完全なinputの取得、GEXとcompanion assayの区別、各streamのread roleの判定です。

UniScFlowでは、このraw dataからの再解析layerを自動化します。

- **accessionに基づく対象選択:** PRJNA/GSEとsample/run関係を解決し、GSM/runの完全一致filterに対応します。元metadataと選択後metadataを両方保存します。
- **複数sourceからの回収:** submitted BAMの有無を確認し、SRA取得と、欠落runに対するENA FASTQ・NCBI source-file fallbackを使います。integrityと選択runのcoverageを検証します。
- **証拠に基づくread解釈:** sample-level GEO protocolをread長、stream間の関係、barcode/UMI信号と照合します。chemistry定義とwhitelistを渡すと10x推定を補強できます。
- **sample単位のrouting:** 明示的GEXをmappingし、明示的non-GEX companionを除外し、未解決の割り当てはreview対象として記録します。全scopeが解決したmixed-platform projectではGSM別に異なる経路を実行できます。
- **source FASTQを書き換えないmapper準備:** canonical linkとsource-to-role対応を作成します。download済みFASTQはmapper preparation時にrename・書換えしません。archiveからの取得時にはadapterがlocal名を付ける場合があります。
- **再現可能なmappingと再開:** STARベースのcommandを生成し、target outputを検証します。完了outputはscopeと実行contextの記録が一致する場合だけ再利用します。localとDockerに対応し、HTML summaryは任意で生成できます。

## Install

repositoryをcloneします。

```bash
git clone https://github.com/KazutoTsukita/scflow-platform-workflow.git uniscflow
cd uniscflow
```

再現可能な利用のため、full Git commit（`git rev-parse HEAD`）と対応environmentまたはcontainer digestを記録してください。release tagは凍結版、`main` はその後の更新を含む版を選択します。

### local commandとして使う

conda環境を作ります。

```bash
conda env create -f environment.yml
conda activate uniscflow
python3 -m pip install -e .
```

開発時にはeditable installが便利です。通常の `python3 -m pip install .` にも対応しており、platform profile、runtime helper script、example configは環境内の `share/uniscflow` にinstallされます。どちらのpip commandも `environment.yml` から作成した環境内で実行してください。pipだけでinstallされるのはPython wrapperと同梱dataであり、STAR、SRA Toolkit、R、GNU Parallel、圧縮toolなどのworkflow実行fileは含まれません。

相対指定したinput/output pathは現在のworking directoryを基準に解決されます。bundled runtime helperとplatform profileはcheckoutまたはinstall済みpackageから自動で解決されます。

以下のcommandで使えます。

```bash
uniscflow --help
```

`environment.yml` で指定されるconda環境には、STAR/STARsolo (`2.7.10b`)、Subreadの `featureCounts`、`samtools`、`wget`、R package、GNU `parallel`、SRA Toolkit、`pigz` などが入ります。

local実行にはBash 4以降とutil-linuxの `flock` commandも必要です。現在のLinux systemには通常ありますが、`environment.yml` からはinstallされません。古いmacOSのsystem Bashはversionが不足し、macOSはdefaultで `flock` を含みません。その場合は、両方を含むDocker imageを使うか、互換host toolを別途installしてください。

### Docker imageとして使う

現在のcheckoutのcodeとtutorialを実行する場合は、同じclean checkoutからimageをbuildしてください。workflow本体、helper script、platform profileを同じsource revisionに揃えられます。

```bash
git clone https://github.com/KazutoTsukita/scflow-platform-workflow.git uniscflow
cd uniscflow
bash docker/build_current_image.sh
```

helperはcleanなGit checkoutだけを受け付けます。Linux/AMD64向けの `uniscflow:latest`、`uniscflow:<version>`、`uniscflow:<7文字commit>` をbuildし、full Git commit、version、build時刻をOCI labelと `/workflow/UNISCFLOW_IMAGE_PROVENANCE` に記録します。CLI、offline plan、25 profile、同梱tool、non-root runtimeもsmoke testします。version文字列だけでは `main` の更新を区別できない場合があるため、full source commitを確認してください。

```bash
docker run --rm uniscflow:latest --version
docker image inspect uniscflow:latest \
  --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
git rev-parse HEAD
```

上の最後の2 commandは同じfull commitを返す必要があります。helperによるprovenance検証とsmoke testを省く最小のlocal buildは以下です。

```bash
docker build --pull --platform linux/amd64 -t uniscflow:latest .
```

source provenanceが必要な場合は `docker/build_current_image.sh` を推奨します。source rebuildはcheckout revisionを固定しますが、upstream base imageやpackage repositoryは変化し得るためbit単位の再現buildではありません。

release workflowのtagged image公開先は `ghcr.io/kazutotsukita/scflow-platform-workflow` です。対応releaseとdigestの公開を確認した場合だけregistry imageを使い、その版のREADME/tutorialに従ってください。現在のcheckoutには上記のsource buildを使えます。full OCI revisionと、公開imageではregistry digestも記録してください。`latest` は便宜的なpointerにすぎません。

imageにはUniScFlow本体、補助tool、STAR/STARsolo、featureCounts、`samtools`、SRA Toolkit、Salmon、R、GNU Parallel、`pigz` が入り、defaultでは非rootのmicromamba userで動きます。genome referenceとlicense付きCell Ranger resourceは含みません。STARsoloとSTAR + featureCounts mappingには、選択したFASTA/GTFから作ったSTAR indexと対応GTFが必須です。Salmonとoptional Cell Ranger targetは、それぞれのindex/referenceを使います。Cell Ranger chemistry定義とbarcode fileはoptionalですが、公共10x dataでは強く推奨します。Cell Ranger mappingを使う場合だけ、Cell Ranger入りの別containerを指定します。

## Quick start

### `--mode all` で最後まで実行する

`--mode all` でmetadata取得、raw-input回収・検証、platform/read構造推定、mapper準備、mappingを実行します。HTML reportには `--write-web-summary` を付けます。以下は10xの例です。小規模な初回実行には、referenceとresource準備も具体的に説明した軽量tutorialを使ってください。

```bash
uniscflow --mode all \
  --ids PRJNA804521 \
  --platform auto \
  --sample-alias GSM5876624,GSM5876625 \
  --filereport-dir work/filereport \
  --download-script-outputdir work/download_script \
  --temporary-sra-download-dir work/sra_tmp \
  --final-file-dir work/raw \
  --mapper-output-dir work/mapper \
  --star-index /path/to/star_index \
  --genes-gtf /path/to/genes.gtf \
  --threads 16 \
  --run-mapper-parallel 1 \
  --max-workers 6 \
  --parallel 6 \
  --cellranger-chemistry-defs /path/to/cellranger/chemistry_defs.json \
  --cellranger-barcodes-dir /path/to/cellranger/barcodes \
  --min-barcode-match-rate 0.8 \
  --inference-report-tsv work/inference_report.tsv \
  --write-web-summary
```

`--sample-alias` を省略すると、そのPRJNA内のすべてのsampleが初期scopeに入ります。その後のsample-level modality/platform routingによりmapper scopeが絞られることがありますが、選択された全sampleの判定は監査記録に残ります。10xでは、Cell Rangerのchemistry定義とbarcode fileを渡すことで、推定chemistryに対応するSTARsolo whitelistを自動選択できます。10x以外では、これらのoptionは必要な場面以外では無視されます。

`--write-web-summary` を付けると、STARsolo sampleではdroplet/UMI用のweb summaryを作り、Smart-seq2のSTAR + featureCounts sampleではfull-length用のweb summaryを作ります。

GEX inference前には、mixed-assay projectのsample-level modality gateを適用します。明示的GEXが1つ以上あり、non-GEXまたはmodality未解決のsampleもある場合、明示的GEXだけをread-structure inferenceとmapper準備へ渡します。明示的non-GEXは除外し、未解決companionはnon-GEXと決めつけずreview対象として保留します。各判定と根拠は `sample_modality_assignment.tsv` に記録します。GEXのmapping成功は、未解決companionも解決したことを意味しません。明示的GEX scopeを確立できない場合、non-GEXだけのscopeはunsupported stop、positiveなmixed-assay evidenceを伴う未解決scopeはreviewのため停止します。

GEX libraryとATAC libraryを明確に分離できる10x Multiome depositでは、GEX成分を通常の10x STARsolo routeでmappingし、ATAC成分をmapper scopeから除外します。GEX sampleには通常の10x chemistry/read-role validationを適用し、ARC-v1と推定された場合は対応するbarcode geometryとwhitelistを使います。明示的GEX scopeを確立できない未解決modalityはreviewのため停止します。

STARsolo GEX reportでは、これとは別にsample multiplexing companionについて高特異度のmetadata-only auditも行います。GEX GSMとは別に `HTO` または `CMO` と明示されたGSMがある場合だけ確定evidenceとし、`cell hashing`、`sample hashing`、`HTO`、`CMO`、`CellPlex` などのsample-level記載はあるもののcompanionを復元できない場合は疑いに留めます。一般的な `multiplexed sequencing`、pooling、CITE-seq、ADT、Feature Barcode、Cell Ranger multi、Illumina indexの記載だけではwarningを出しません。multiplex判定自体はwarning-onlyであり、mapping成功判定、halt logic、exit statusを変更しません。UniScFlowは対象となるGEX成分をmappingしますがHTO/CMOによるsample demultiplexingは行わないため、flagされたmatrixはpooled sampleを含む可能性があり、projectのmanual reviewが必要です。

### `--mode all` が行うこと

```text
1. project metadataを解決し、指定sample/run filterを適用する。
2. currentな検証済みinputを再利用し、欠落runを対応sourceから取得する。
3. input integrityと選択run全体のcoverageを検証する。
4. sample-level assay metadataとraw-read evidenceを照合する。
5. 対象sampleのplatform、chemistry、logical read roleを推定する。
6. canonical FASTQ link、または適格なtag付きBAM inputを準備する。
7. 適切なSTARsoloまたはSTAR + featureCounts commandを生成・実行する。
8. mapper outputを検証し、scopeに結び付いた完了receiptを記録する。
9. inference記録、manifest、log、任意のHTML summaryを出力する。
```

BD Rhapsody、Parse/Evercode、SPLiT-seq、sci-RNA-seq、CEL-seq/MARS-seq、inDropsなど、dataset固有のbarcode manifestやwell manifestが必要なplatformを検出した場合、`--mode all` はdownloadと推定後に意図的に停止します。これは、それらしいが誤ったmapping結果を作らないための安全機構です。

mapper完了判定はdefaultでstrictです。mapper scopeに割り当てられた全SRRがintegrity check済みFASTQ/BAM inputに対応し、生成された全sample mapperがexit 0に加えて、空でない構造的に有効なendpoint matrixを生成した場合だけmapped endpointになります。欠落run、失敗sample、またはoutput欠損があればnon-zeroを返し、manifestに記録します。`--allow-partial-success` は一部sampleのoutput欠損を明示的に許容する場合だけ使い、project全体の完全な結果が必要な場合には使わないでください。rerunでは、current scopeをすでに満たす検証済みinputをdownload前に再利用します。意図した `documented_halt`、`non_target_stop`、`unsupported_stop` は正常なterminal endpointであり、matrixを作らずexit 0になることがあります。`needs_review` は未解決状態でありnon-zeroを返します。

### 終了結果の読み方

**documented halt (recognized stop)** は内部endpointの `documented_halt` に対応します。Recognized stopとunsupported stopは別カテゴリーで、両方を **actionable stop** に含めます。

| 結果 | 確認方法 | 意味 |
| --- | --- | --- |
| Mapped | `mapper_run_manifest.tsv` が各mapper対象sampleを `ok` または `reused` と記録し、target固有のoutput検証に合格し、各mapper-input directoryにcurrentな `.uniscflow_mapping_complete.json` がある | そのmapper scopeのmappingが完了しています。`sample_modality_assignment.tsv` と、存在する場合は `sample_platform_routing.tsv` を確認し、mapped sample、除外companion、その他のrouteを区別します。 |
| Documented halt (recognized stop) | current halt markerに `halt_type=manual_preprocessing_required` とprofile固有guidanceがある。必要に応じて `halt_summary.html` も作られる | platformは判定済みですが、外部manifest、barcode、vendor preprocessingが必要です。matrix完了を意味しません。 |
| Unsupported stop | halt markerの `halt_type` が `unsupported_platform` または `non_target_data` | assayは認識済みですが検証済みUniScFlow mapping routeがないか、通常のbulk RNA-seqなどsc-/snRNA-seq対象外のscopeです。matrix完了を意味しません。 |
| Needs review | inferenceまたはroutingの進行を妨げる `needs_review` が記録され、commandがnon-zeroを返す | 要求されたrouteのevidenceが不完全または矛盾しています。再実行前に根拠を確認してください。明示的GEX scopeから保留された未解決companionは、上記のとおり別途記録されます。 |
| Temporary metadata service | exit statusが `75` | ENA/GEOの一時障害または繰り返し壊れたtransport responseです。生物学的endpointは未判定なので、後で同じcommandを再実行します。 |
| その他の失敗 | project logと `mapper_run_manifest.tsv` を確認する | input、integrity、inference、mapperの報告errorを解決します。最終matrixに `--allow-partial-success` を使わないでください。 |

STARsolo droplet経路のendpoint検証では、空でない構造的に有効なraw GeneまたはGeneFull matrixと、対応log・summaryを確認します。raw matrixの列はbarcodeであり、すべてがcalled cellとは限りません。mapping成功は非空のfiltered-cell matrixを保証せず、後続のcell-quality評価の代わりにもなりません。Smart-seq2の行列列は、推定cell granularityまたは指定grouping mapに従います。

保存されたplatform reportとhalt markerは、選択sample alias、SRR accession、filtered filereportのSHA-256、およびinput FASTQ/BAMの相対path・size・mtime/ctimeから作るfingerprintに結び付けられます。別scopeや、同じpathで後から置換されたinputの記録は再利用されません。mapping receiptはtarget、reference context、command checksumにも結び付けられ、再利用前にoutputを再検証します。また、状態を変更するcommandではPRJNA単位のfilesystem lockを保持します。同じprojectを誤って二重起動した場合は、manifestやmapper outputを相互上書きせず直ちに停止します。

1 commandに複数のPRJNA/GSE IDを渡した場合も、platform inferenceとmapper preparationはBioProjectごとに独立したruntime configで実行されます。先のprojectで推定されたplatformが次のprojectへ引き継がれることはありません。

### 補助mode

以下のmodeは、環境確認、debug、途中からの再開、中間結果の確認に使う補助modeです。UniScFlowの基本workflowは `--mode all` です。

```bash
uniscflow --mode check ...
uniscflow --mode download ...
uniscflow --mode infer-platform ...
uniscflow --mode mapping ...
uniscflow --mode report ...
```

### Dockerで最後まで実行する

以下ではlocal buildした `uniscflow:latest` を使います。保存解析ではbuild helperが表示したcommit-specific tagへ置き換え、公開済みcommit tagは移動させないでください。work directoryはread-write、対応するSTAR index/GTFはread-onlyでmountします。

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v /path/to/filereport:/data/filereport \
  -v /path/to/download_script:/data/download_script \
  -v /path/to/sra_tmp:/data/sra_tmp \
  -v /path/to/raw_fastq:/data/raw_fastq \
  -v /path/to/mapper:/data/mapper \
  -v /path/to/reference:/ref:ro \
  uniscflow:latest \
  --mode all \
  --ids PRJNA804521 \
  --platform auto \
  --sample-alias GSM5876625 \
  --filereport-dir /data/filereport \
  --download-script-outputdir /data/download_script \
  --temporary-sra-download-dir /data/sra_tmp \
  --final-file-dir /data/raw_fastq \
  --mapper-output-dir /data/mapper \
  --star-index /ref/star_index \
  --genes-gtf /ref/genes.gtf \
  --threads 16
```

imageはdefaultでも非rootのmicromamba userで動作します。hostのUID/GIDを明示すると、mountしたwork directoryへ書かれるfileもhost user所有になります。`HOME=/tmp` は、numeric host userで動くtoolがuser設定やcacheを作るための書き込み可能なhomeを提供します。SRA、展開後FASTQ、STAR一時file、mapper outputが同時に存在し得るため、十分な空き容量を確保してください。

proxyが必要な環境だけ `--ftp-proxy` を追加します。通常は指定不要です。

```bash
--ftp-proxy proxy.example.ac.jp:8080/
```

## mode

`all` が主役のmodeです。その他のmodeは、debug、監査、途中からの再開、reference作成などのためにworkflowの一部を取り出した補助modeです。

| mode | 役割 | 内容 |
| --- | --- | --- |
| `all` | 基本workflow | inputを取得・検証し、platform/read roleを推定して、適格なmapping経路を実行するかactionable stopを報告します。HTML summaryは任意です。 |
| `check` | 補助 | full runの前に、必要なcommand、package、script、mapping設定を確認します。 |
| `plan` | 補助 | workflow stageを起動せず、実行予定のcommandを表示します。設定directoryは作成されることがあり、GSE inputではGEO metadataの解決とcache書き込みが行われることがあります。 |
| `download` | 補助 | metadata取得、複数sourceからのinput回収、integrity確認、platform/read推定を行い、mappingは実行しません。mapper準備前のinput・推定記録の確認に使います。 |
| `infer-platform` | 補助 | 既存metadataとraw inputからplatform推定を再実行します。 |
| `infer-reads` | 補助 | 既存の `SRR*_[0-9].fastq.gz` directoryからread構造だけを推定します。 |
| `prepare` | 補助 | 既存FASTQまたは適格なtag付きBAMからmapper input、script、manifestを生成し、mappingは実行しません。 |
| `mapping` | 補助 | 選択routeに従って既存inputからmapper scriptを準備・実行します。raw dataのdownloadは行いません。 |
| `report` | 補助 | 既存のSTARsoloまたはSTAR + featureCounts mapper outputから `web_summary.html` を生成します。 |
| `validate` | 補助 | metadata、生成script、raw-input directory、選択run coverageを含むinput準備状態を確認します。 |
| `build-star-index` | utility | FASTA/GTFからSTAR genome indexを作成します。 |

shortcutも使えます。

```bash
--modecheck
--modedownload
--modemapping
--modeall
```

## よく使うoption

この表は通常利用で必要なoptionをまとめたものです。完全かつ最新のCLI referenceは `uniscflow --help` で確認してください。

| option | 内容 |
| --- | --- |
| `--ids` | PRJNAまたはGSE project ID。BioProject IDでは`PRJNA` prefixはあってもなくてもよいです。GSE IDはGEO SOFTの明示的なBioProject relation fieldだけからPRJNAへ解決されます。 |
| `--sample-alias` | GSMでfilterします。省略するとproject全体が初期scopeに入り、その後modality/platform routingでmapper scopeが絞られる場合があります。複数指定はカンマ区切りです。 |
| `--sample-map-tsv` | plate/full-length project用のgrouping tableです。GSM列は `gsm_accession`, `gsm`, `sample_alias`, `geo_accession`、生物学的sample列は `sample_id`, `biological_sample`, `sample_name`, `group_id`, `sample` のいずれかです。任意のcell/well列は `cell_id`, `well_id`, `well`, `cell`, `library_id`、condition列は `condition`, `treatment`, `group` を使えます。 |
| `--run-accession` | SRRやERRなどのrun accessionを完全一致でfilterします。 |
| `--filter` | ENA metadata列を `key=value` でfilterします。複数回指定できます。 |
| `--filereport-dir` | ENA metadataの保存先です。 |
| `--download-script-outputdir` | 生成したSRA download scriptの保存先です。 |
| `--temporary-sra-download-dir` | SRA fileと中間FASTQの一時directoryです。 |
| `--final-file-dir` | download FASTQとsubmitted BAMのraw-input保存先です。mapper preparationはsource FASTQを変更せず、別途 `--mapper-output-dir` 以下にcanonical linkを作ります。 |
| `--max-workers` | `fasterq-dump` の並列worker数です。 |
| `--parallel` | download、変換、圧縮で使う並列度です。 |
| `--ftp-proxy` | FTP proxyです。通常は不要で、大学などproxyが必要な環境だけ指定します。 |
| `--platform` | `auto` でmetadataとraw-read evidenceを照合します。`10x`, `dropseq`, `smartseq2` などを明示しても、その要求とevidenceの整合性を確認します。 |
| `--force-platform` | metadata/FASTQ矛盾に対するexpert向けoverrideです。non-Flex haltを解除できる場合がありますが、不足manifestやvendor resourceを補うものではなく、正しいmappingの証拠にもなりません。専用の `10x_flex` 保護は解除できません。 |
| `--generic-cell-barcode-read`, `--generic-cell-barcode-start`, `--generic-cell-barcode-length`, `--generic-umi-read`, `--generic-umi-start`, `--generic-umi-length`, `--generic-cdna-read` | `--platform generic_droplet_umi` 専用の完全な明示geometryです。7項目すべてが必須です。cell barcodeとUMIは同じlogical read上の重ならない区間、cDNAはもう一方のlogical readでなければなりません。UniScFlowはこれらを自動推測しません。 |
| `--resolve-bam` | defaultで有効です。ENA/GEO上でraw split FASTQではなくsubmitted BAM/alignment fileが登録されている場合に、BAMを直接downloadし、barcode/UMI SAM tagを検査して、完全なraw `CR/CY/UR/UY` tag（または由来を確認した旧Cell Rangerの`CR/CQ/UR/UQ` tag）がある場合だけSTARsolo-from-BAM rescue scriptを作ります。 |
| `--no-resolve-bam` | 初期のsubmitted-BAM rescue検出を無効化します。SRA取得と欠落runのsource fallbackは引き続き使えます。 |
| `--bam-integrity-check` | submitted BAM rescue fileを受理する前の検証方法です。`full` は `samtools view -c` で全走査し、CRC/inflate errorを検出できます。`quickcheck` は速いですが弱いです。defaultは `full` です。`samtools` がなければfail-closedで停止し、ENA byte数とMD5も値があれば照合します。 |
| `--bam-integrity-retries` | BAM integrity検証に失敗したsubmitted BAMを追加でdownloadし直す回数です。defaultは `2` です。 |
| `--fastq-integrity-check` | ENA `fastq_ftp` やNCBI SDL source fileから直接downloadしたFASTQ.gz fallback fileを受理する前の検証方法です。`gzip` は圧縮file全体を読み、gzip CRC/truncation errorを検出し、ENA byte数とMD5も値があれば照合します。`none` でgzip scanを無効化できます。defaultは `gzip` です。 |
| `--fastq-integrity-retries` | FASTQ.gz integrity検証に失敗したdirect FASTQを追加でdownloadし直す回数です。defaultは `2` です。 |
| `--profiles-dir` | platform profile JSON fileのdirectoryです。通常はbundledの `profiles/platforms` を使うため指定不要です。 |
| `--geo-soft-max-samples` | 初期GEO metadata summaryに使う先頭GSM数です。defaultは `3` で、選択sample全体のevidence確認の上限ではありません。 |
| `--geo-soft-dir` | GEO SOFT cache保存先です。省略時は `<filereport-dir>/geo_soft` を使います。 |
| `--target` | `--mode prepare` のtargetです。実行・検証可能なmapper targetは `starsolo`, `star_featurecounts`, `cellranger`, `salmon` で、`auto` はprofile defaultを選びます。`manual_review` はmapper commandではなくreview artifactを書く非実行preparation targetです。明示したtargetは `--mapping-engine` より優先され、legacy Cell Ranger engineはDocker routeを選びます。外部vendor/third-party ecosystem候補はprofileの `external_targets` に分離され、placeholder commandへ黙ってfallbackしません。 |
| `--mapper-output-dir` | platform-awareなmapping scriptやmanifestの出力先です。 |
| `--star-index`, `--salmon-index` | STARsolo、STAR + featureCounts、optional Salmon script生成時に使うindexです。 |
| `--genome-fasta`, `--genes-gtf`, `--sjdb-overhang` | `--mode build-star-index` で使う入力です。`--sjdb-overhang` のdefaultは `99` です。`--genes-gtf` はSmart-seq2のSTAR + featureCounts mappingでも使います。 |
| `--starsolo-whitelist` | 生成されるSTARsolo scriptに渡すbarcode whitelistを手動上書きするoptionです。10xでは `--cellranger-chemistry-defs` と `--cellranger-barcodes-dir` があれば、推定chemistryから自動選択できます。 |
| `--mapping-engine` | `--mode mapping` のrunner familyです。defaultの `starsolo` は、STARsoloとSTAR + featureCountsを含む生成済みscriptを実行します。`cellranger` と `--target auto` の組合せはoptional legacy Docker routeです。明示的な `--target cellranger` はmapper scriptを生成し、hostから実行できる `cellranger` executableを必要とします。 |
| `--run-mapper-parallel` | 生成済みmapper scriptを同時にいくつ実行するかです。 |
| `--allow-partial-success` | 一部sample mapperが失敗しても、1件以上成功すればproject commandを成功扱いする明示的opt-inです。defaultは全sample成功を要求します。 |
| `--write-web-summary` | mapping後、STARsoloまたはSTAR + featureCounts outputがあるsampleに `web_summary.html` を生成します。 |
| `--report-name` | 生成されるHTML reportのfile名です。defaultは `web_summary.html` です。 |
| `--no-logo` | UniScFlow起動時のlogoを表示しません。logoはterminal表示のみで、project logには書き込まれません。 |
| `--version` | install済みUniScFlow versionを表示します。 |
| `--index1`, `--index2`, `--read1`, `--read2` | FASTQ suffixとI1/I2/R1/R2の対応を手動指定します。ないindex readは `NULL` にします。 |
| `--cellranger-chemistry-defs` | 10x chemistry-aware推定に使うCell Rangerの `chemistry_defs.json` です。 |
| `--cellranger-barcodes-dir` | 10x chemistry-aware推定に使うbarcode whitelist directoryです。 |
| `--min-barcode-match-rate` | chemistry/read-role推定で要求するbarcode一致率です。defaultは `0.5` で、例の `0.8` は意図的に厳しい指定です。scope、assay、geometryの検証も別途必要です。 |
| `--cellranger-chemistry` | 推定対象のCell Ranger chemistryを絞ります。複数回指定できます。 |
| `--inference-report-tsv` | project-level 10x read-inference pathが出力する場合に、FASTQごとの詳細をTSVへ追記します。全platform共通のinference logではなく、deferred non-10x routeはsample別assignment fileを使います。 |
| `--cellranger-container` | Cell Rangerが入っている既存Docker container名です。 |
| `--transcriptome` | Cell Ranger container内のreference transcriptome pathです。 |
| `--localcores`, `--localmem` | Cell Rangerに渡す計算資源設定です。 |
| `--no-bam` | mapper BAM outputを残しません。defaultの `auto` policyでは、STARsolo targetはBAM outputを避け、STAR + featureCountsはcount生成に必要な一時BAMだけを書き、count生成後に削除します。 |
| `--with-bam` | 対応targetでmapper BAM outputを保持します。通常はcount-onlyになるSTARsolo targetや、通常は削除されるSTAR + featureCountsの一時BAMも保持します。 |
| `--cellranger-include-introns` | legacy `cellranger count` に明示的な `true` または `false` を渡します。省略時はCell Ranger versionごとのdefaultを使います。 |
| `--config` | TOML config pathです。同梱exampleはすでにdefaultで、checkout内では `./config/example.toml` をtemplateとして使えます。 |
| `--dry-run` | workflow stageを起動せず、予定actionを表示します。設定directoryは作成されることがあり、GSE解決ではmetadataをdownloadしてcacheすることがあります。 |
| `--help` | 完全かつ最新のoption referenceを表示します。 |

異なるkeyのmetadata filterは **AND** で結合され、1つのfilter内のカンマ区切り値は **OR** で結合されます。同じkeyを繰り返した場合は最後の値だけが残ります。GSM/sample aliasやSRR run accessionなどのaccession列は大文字小文字を無視した完全一致で判定するため、`GSM1`が`GSM10`を選ぶことはありません。`sample_title`などの自由文列は、大文字小文字を無視した固定文字列の部分一致です。存在しない列や不正なfilterはdownload前に停止します。

## 対応platform

UniScFlowでは、25個の名前付きplatform routing profileを定義しています。これとは別に、platform profileではなくユーザー指定geometryとして扱うconfigurable generic droplet-UMI routeがあります。以下では、自動mapping、resource依存のdocumented halt (recognized stop)、認識済みunsupported assayを分けて示します。

| level | platform/profile | 現在の対応 |
| --- | --- | --- |
| 自動mapping | 名前付きplatform: `10x`, `dropseq`, `seqwell`, `smartseq2`; optional user-configured route: 完全な明示geometryを伴う `generic_droplet_umi` | mapper script/manifestを生成して実行できます。10x、Drop-seq、Seq-Well、明示設定されたgeneric droplet-UMI routeはSTARsoloです。defaultのSmart-seq2 STAR + featureCounts mappingはcurrent scopeに一致するgranularity auditを行い、GSM-as-cell inputはSTAR + featureCounts、review済みsample-map groupまたは自動判定されたrun-as-cell inputはSTARsolo SmartSeqを使い、ambiguous layoutは停止します。明示的に選んだSalmon targetはこのgranularity routingを通らないため、GSM/cell/sampleの意味を手動確認した場合だけ使用します。10xのCell Ranger script生成もoptional legacy targetとして残します。 |
| Documented halt (recognized stop) | `hive_clx`, `bdrhapsody`, `bdrhapsody_targeted_panel`, `dnbelab_c4`, `parse`, `splitseq`, `scirnaseq`, `celseq2`, `marsseq`, `indrop`, `scrbseq`, `smartseq3`, `microwellseq`, `singleron_gexscope`, `seekone`, `mobidrop_mobicube`, `fluidigm_c1`, `icell8`, `ramda_seq`, `quartz_seq`, `pipseq` | downloadとplatform推定までは行い、実験固有のbarcode geometry、well/sample map、vendor preprocessing、manifestが必要なためhalt markerを書いて停止します。各markerにはprofile固有の停止要因、必要な追加入力、具体的な次の手順、推奨workflow、再開条件が入り、CLIとhalt用web summaryにも同じ内容を表示します。標準CLIでは、halt markerがある間は後続の `prepare` や `mapping` もそのsampleをskipします。詳細は[documented halt (recognized stop) profiles: what to do next](docs/recognized_halt_next_steps.md)を参照してください。 |
| Unsupported stop | 通常のbulk RNA-seq、targeted transcriptomics、Stereo-seq/spatial transcriptomics、ddSEQ、ATAC-onlyなど、認識済みのnon-targetまたはunsupported assay | exact selected scopeとevidenceを記録してhalt markerを書き、mapper outputを作らず停止します。読者向けにはUnsupported stopとし、内部記録では `non_target_stop`（`halt_type=non_target_data`）と `unsupported_stop`（`halt_type=unsupported_platform`）を区別します。 |

PIP-seqは `pipseq` として認識しますが、PIP-seq chemistry固有のbarcode処理をgeneric STARsolo設定で代用しません。正確なchemistry/barcode resourceを使うPIPseeker互換workflowで処理し、検証済みmatrixを下流解析へ渡してください。

10x Flexは通常の10x STARsoloではなく、専用のprotected halt routeです。対応probe set、sample/probe-barcode設定、`cellranger multi` または別の検証済みprobe-aware workflowが必要です。`--force-platform` でもこの保護は解除できません。

Stereo-seqは、名前付きprofileではない `spatial_transcriptomics` routing labelとして認識します。`--platform auto` では、選択された各GSMにsample-levelのStereo-seq/STOmics assay evidenceと、対応する `SAW`、`Stereo-seq Analysis Workflow`、または `Stereopy` のprocessing/output evidenceが必要です。versionのない大文字 `SAW` は受理しますが、assayだけ、processingだけ、通常の小文字 `saw`、Seriesだけ、一部GSMだけのevidenceでは不足です。確定routeは `unsupported_stop`、`halt_type=unsupported_platform` で終わり、mapper outputは生成しません。

platform familyが混在するBioProjectでは、選択された各GSMが完全、重複なしのendpointを1つだけ持ち、`needs_review` が残らない場合にsample-level routeが成功します。mapper commandは `automatic_mapping` sampleだけに生成し、`documented_halt`、`non_target_stop`、`unsupported_stop` sampleは明示的terminal endpointとして残します。重複、不完全、overlap、needs-review planは停止し、より狭いselectionまたは追加evidenceが必要です。全sampleが10xでchemistryやread-role layoutだけがGSMごとに異なる場合は、mapper preparation時のsample-level 10x inferenceで自動正規化できます。

Plate-basedのlibrary構築法名はbulk RNA-seq depositにも現れます。このためUniScFlowは、すべてのplate-family profileについてtechnology候補を同定した後に、高特異度のnon-target guardを適用します。`bulk RNA-seq`、`bulk, RNA-seq`、`bulk 3'-end RNA-seq`、`bulk 3-prime RNA-seq`、`bulk transcriptome/transcriptomics` など、assayを明示する記述があれば `non_target_bulk_rna` にroutingします。sample-levelの明示記述は単独で十分ですが、Series-levelの記述ではgene-by-sample matrix、sample barcode、sampleごとのdemultiplexingなど、barcodeがcellではなくsampleを表す根拠も必要です。platform名を解決できなくても、FASTQがplate/full-length-likeと判定された場合には同じpatternを追加評価します。このFASTQ-gated routeには、transcriptomic RNA-seq metadata、supported barcode/UMI geometryがないこと、明示的bulk RNA-seq evidence、sample-levelおよびSeries-levelのsingle-cell evidenceがないことをすべて要求し、推定済みのnon-plate platformを上書きしません。推定されたplate technologyがある場合はreportとhalt markerの `technology_candidate` に監査情報として残しますが、mapper commandやcount matrixは生成しません。単独の `bulk`、homogenized tissue、FACS sorting、low input、またはplatform比較への言及だけでは、このplate-family共通vetoは発動しません。

DNBelab C4/C Series系の3' scRNA-seqでは、DNBelab C4、DNBSEQ/DIPSEQ、PISA、DNBC4Toolsなどのmetadata evidenceと、DNBelabらしいbarcode/UMI readというFASTQ evidenceが揃う場合に `dnbelab_c4` profileを自動選択します。ただしこれはplatform-aware routingであり、defaultの自動mapping対象ではありません。DNBelab/PISA互換のbarcode補正やpreprocessingが必要になり得るため、`--mode all` はdownload後にhalt markerを書いて停止し、STARsolo commandを直接生成しません。C Seriesなどの製品family名だけでは強く決め打ちしません。

Honeycomb Biotechnologies HIVE CLX depositは、HIVE CLXまたはBeeNetの根拠がある場合に `hive_clx` として記録します。deposit内でcommercial Seq-Well platformと記載されていても、通常のSeq-Well STARsolo profileには送りません。BeeNet互換のvendor preprocessingとbarcode geometryが必要なため、`--mode all` は理由を記録したmanual-preprocessing haltで停止します。

BD Rhapsody targeted-expression depositは、選択したすべてのGSMがBD Rhapsody Immune Response Panel、targeted workflow、明示的なgene/transcript数を個別に示し、WTA evidenceがない場合に限り `bdrhapsody_targeted_panel` として記録します。このsubtypeも理由付きhaltとなり、正確なpanel/custom-primer target、対応するBD barcode/reference resource、使用時のsample-tag/AbSeq定義、BD Rhapsody Targeted Analysis Pipelineが必要です。V(D)J enrichmentを含む一般的なtargeted-PCR記載だけではこのrouteを選びません。

Singleron GEXSCOPE/SCOPE-chip datasetでは、Singleron、GEXSCOPE、SCOPE-chip、GEXSCOPE Single-Cell RNA Library Kitなどのmetadata語を検出します。これはvendor-specific droplet UMI dataとして扱います。公共FASTQは150 bp paired-endのように見えることがありますが、barcode/UMI抽出やpoly-A/adapter処理がSingleron preprocessing workflowに含まれるため、UniScFlowはplatformを記録したうえでdownload後に停止します。validatedなSingleron preprocessing pathを追加するまでは、自動mappingには進みません。

HIVE CLX/BeeNet、BD Rhapsody、DNBelab C4/C Series、Parse/Evercode、SPLiT-seq、sci-RNA-seq、CEL-seq/MARS-seq、inDrops、SCRB-seq、Microwell-seq、SeekOne、Smart-seq3、MobiDrop/MobiCube、Fluidigm C1、ICELL8、RAMDA-seq、Quartz-seqなど、実験固有のbarcode/manifestやvendor preprocessingが必要なplatformでは、`--mode all` は通常downloadとplatform推定の後で停止します。project FASTQ directoryに `.uniscflow_halt_after_download.json` を書き、誤った結果を作らないようにmapper command生成とmappingをskipします。active markerがあるscopeは後の `prepare` と `mapping` でもskipされ、resource fileを置くだけでは自動解除されません。`resume_mode=external_workflow` のprofileは指定external pipelineで完了します。`manual_review_then_uniscflow` はmarkerの案内に従い、review後に `--force-platform` と明示的supported mapper targetを使った場合だけ継続できます。低レベルのmapper-input helperをmanual-review target/profileで手動実行し、review artifactとsource-role symlink FASTQを作成できる場合もあります。

通常のUMI-bearing Smart-seq3はdocumented halt (recognized stop)です。`--platform auto` の場合だけ、modified Smart-seq3 protocolが、選択された全GSMについてnon-UMI、one-cell-per-well、complete input、competing assayなしという限定checkをすべて満たせば、Smart-seq2 computational backendを使うことがあります。それ以外ではSmart-seq3 haltを維持します。

各profileの `supported_targets` には、UniScFlowが実行commandを生成しendpointまで検証できるtargetだけを記録します。manual preprocessing後に利用候補となるvendor/third-party pipelineは、実行可能targetと混同しないよう `external_targets` に分離しています。

Smart-seq2以外のnon-droplet platformと推定され、かつ選択されたGSM/sample_aliasが96以上ある場合、UniScFlowはproject logに太字のwarningを出します。`--mode all` ではdownload後にいったん停止し、mappingへは進みません。Smart-seq2はdefaultのSTAR + featureCounts routeが上記のpost-inference granularity auditを行うため、このgeneric stopの対象外です。明示的Smart-seq2 Salmon targetはそのroutingを通らないため、GSM/cell/sampleの意味を手動確認する必要があります。plate系やwell-indexed系の公共datasetでは、1つのGSMが1つのwell/cellを表し、生物学的sampleではないことがあるため、sample-level outputを解釈する前にGEO/SRA metadataを手動確認してください。確認後に続ける場合は、明示的に `uniscflow --mode mapping` を再実行するか、`--sample-map-tsv` を渡してGSM/well outputを生物学的sampleごとにまとめます。

## 前提条件

### 好きなFASTA/GTFからSTAR indexを作成する

STARsoloとSTAR + featureCounts mappingでは、STAR genome indexと対応GTFが必要です。均一に再解析したいreference genome FASTAとgene annotation GTFを選び、最初に一度だけSTAR indexを作成します。Salmonとoptional Cell Ranger targetは、それぞれのindex/referenceを使います。

```bash
uniscflow --mode build-star-index \
  --star-index /path/to/star-index \
  --genome-fasta /path/to/genome.fa \
  --genes-gtf /path/to/genes.gtf \
  --threads 16
```

10x風のreference directoryを使う場合は、例えば以下のようになります。

```bash
uniscflow --mode build-star-index \
  --star-index /path/to/star-mm10-2020-A \
  --genome-fasta /path/to/refdata-gex-mm10-2020-A/fasta/genome.fa \
  --genes-gtf /path/to/refdata-gex-mm10-2020-A/genes/genes.gtf \
  --threads 16
```

UniScFlowはSTAR index作成時に、defaultで `--sjdb-overhang 99` を使います。公共scRNA-seq datasetを広く再解析するための実用的なdefaultです。datasetごとに厳密に最適化したい場合は、STARの推奨通り `read length - 1` を明示的に指定できます。作成後は、`--mode all` に同じpathを渡します。

```bash
--star-index /path/to/star-mm10-2020-A
--genes-gtf /path/to/refdata-gex-mm10-2020-A/genes/genes.gtf
```

`build-star-index` はFASTA/GTFのSHA-256 provenanceを `uniscflow_star_index_manifest.json` に記録します。出力path専用lockを保持し、一時的な同階層directoryで完全なindexを構築・検証してから設定pathへ昇格するため、再構築が失敗しても従来indexは保持されます。STAR-based mapping前には、`--genes-gtf` をこの記録と照合します。STARを直接実行して作った既存indexでは、`genomeParameters.txt` の `sjdbGTFfile` が現在も参照可能ならそのGTFをhash照合します。hard errorにするのは、GTF hashの明示的不一致を確認できた場合だけです。provenanceがない、読めない、古い、または過去のGTF pathへ到達できない場合は、それ自体を不一致とは扱いません。STAR indexと現在のGTFのgene ID overlapを報告し、`star_index_annotation_provenance_unverified` warning付きでmappingを続け、そのwarningをmapper profileとweb summaryにも残します。最終的なproduction解析では、完全な監査性のためUniScFlowでindexを再構築することを推奨します。

### Optionalだが推奨: Cell Ranger chemistry/barcode file

UniScFlowのdefault mappingはSTARベースなので、Cell Ranger本体は必須ではありません。ただし公共10x datasetでは、Cell Rangerの `chemistry_defs.json` とbarcode inclusion listを渡すと、read構造・chemistry推定の精度がかなり上がります。UniScFlowは、sampleしたFASTQ barcode readをchemistry-specific barcode listに照合し、STARsoloに渡すwhitelistも自動選択します。

STARsolo script生成時には、profileのCB/UMI終端まで届くreadの割合をbarcode-read FASTQごとに確認します。sampled readの90%以上なら通常どおり進み、70%以上90%未満ならshort-barcode-read warningを残します。70%未満では1,000、5,000、10,000 readsへsamplingを増やし、それでも不足する場合は、下記の限定的な短縮UMI経路に該当しない限り停止します。必要範囲より長いreadは正常であり、長さだけで別chemistryとは判定しません。aggregate chemistry推定がactionableな場合もFASTQ別whitelist一致率warningを保持し、mapper profileと任意のHTML summaryへ表示します。

**均一に短い10x UMI:** 対応するSC3Pv3、SC3Pv4、`SC5P-R2-v3`、`ARC-v1`定義では、chemistryが `min_length` で明示的に許容する場合に限り、公称長より短い観測UMI長を使えます。完全な16-base barcode、単純なR1 barcode/UMI + R2 cDNA構造、各barcode FASTQで十分なwhitelist一致、sample内の全barcode streamが同じread長であることの全file走査確認が必要です。source FASTQやchemistry定義は変更せず、有効geometryを記録します。`shortened_umi` warningには、UMI塩基数の減少により分子の識別衝突が増える可能性を示します。長さの混在や非対応geometryには適用しません。

10x Flex/Fixed RNA Profiling chemistry（`MFRP-*`および`Flex-v2-*`）は、`10x_flex` routing subtypeとして記録します。FlexはSTARsoloへ流さず、public inputのdownloadと整理後に理由付きで停止します。専用のCell Ranger `multi` workflowには、対応するprobe setとsample/probe-barcode設定が必要です。

UniScFlowが必要とするCell Ranger resourceは以下だけです。

| 必要なresource | Cell Ranger展開後の典型的な場所 | 用途 |
|---|---|---|
| `chemistry_defs.json` | `cellranger-x.y.z/lib/python/cellranger/chemistry_defs.json` | Cell Ranger chemistry名、read role、barcode長、offset、whitelist file名を読むために使います。 |
| `barcodes/` directory | `cellranger-x.y.z/lib/python/cellranger/barcodes/` | `3M-february-2018*.txt.gz`, `737K-august-2016*.txt`, ARC/ATAC/GEX系などの10x barcode inclusion listを含みます。 |

`barcodes/` は1つのwhitelist fileだけではなく、directoryごとcopyしてください。公共depositでは10xの世代がさまざまなので、UniScFlowはchemistry候補ごとに対応するbarcode listへ照合します。

10x Genomics公式download pageからCell Rangerをdownloadし、10xのlicense termsに同意した上でtarballを展開します。最新のversion、file名、download URLは必ずCell Ranger websiteで確認してください。

https://www.10xgenomics.com/support/software/cell-ranger/downloads

```bash
mkdir -p "$HOME/cellranger_download"
cd "$HOME/cellranger_download"

# このblockの前に、10x websiteが表示する現在のsigned URLをexportします。
# 例: export CELLRANGER_URL='<current signed Cell Ranger URL>'
: "${CELLRANGER_URL:?CELLRANGER_URLに10xの現在のsigned URLを設定してください}"

wget -O cellranger-10.0.0.tar.gz "$CELLRANGER_URL"
tar -xzvf cellranger-10.0.0.tar.gz
```

Cell Ranger 10.0.0を展開した場合、UniScFlowに指定する2つの入力は以下です。

```text
cellranger-10.0.0/lib/python/cellranger/chemistry_defs.json
cellranger-10.0.0/lib/python/cellranger/barcodes/
```

別のCell Ranger versionを使う場合は、`10.0.0` をdownloadしたversionに置き換えてください。

```text
cellranger-x.y.z/lib/python/cellranger/chemistry_defs.json
cellranger-x.y.z/lib/python/cellranger/barcodes/
```

UniScFlowで使いやすいように、必要なfileを小さな固定directoryにcopyしておくと便利です。

```bash
CELLRANGER_HOME="$HOME/cellranger_download/cellranger-10.0.0"

mkdir -p "$HOME/cellranger_v10_for_scflow/barcodes"
cp "$CELLRANGER_HOME/lib/python/cellranger/chemistry_defs.json" \
  "$HOME/cellranger_v10_for_scflow/chemistry_defs.json"
cp -R "$CELLRANGER_HOME/lib/python/cellranger/barcodes/." \
  "$HOME/cellranger_v10_for_scflow/barcodes/"

CELLRANGER_CHEMISTRY_DEFS="$HOME/cellranger_v10_for_scflow/chemistry_defs.json"
CELLRANGER_BARCODES_DIR="$HOME/cellranger_v10_for_scflow/barcodes"

ls "$CELLRANGER_CHEMISTRY_DEFS"
ls "$CELLRANGER_BARCODES_DIR"
```

Cell Ranger packageの内部layoutが違う場合は、先に以下で場所を探してください。

```bash
find "$CELLRANGER_HOME" -name chemistry_defs.json
find "$CELLRANGER_HOME" -type d -name barcodes
```

準備後のresource directoryは以下のようになります。

```text
$HOME/cellranger_v10_for_scflow/
  chemistry_defs.json
  barcodes/
    3M-february-2018.txt.gz
    3M-february-2018_TRU.txt.gz
    737K-august-2016.txt
    ...
```

`--mode all` では以下のように渡します。

```bash
--cellranger-chemistry-defs "$CELLRANGER_CHEMISTRY_DEFS"
--cellranger-barcodes-dir "$CELLRANGER_BARCODES_DIR"
--min-barcode-match-rate 0.8
```

これらは推定とwhitelist選択のためのresourceとして使われます。UniScFlowのdefault STAR-based workflowではCell Ranger本体は実行しません。また、UniScFlowはCell Ranger本体や10x barcode fileを再配布しません。

## 出力ファイル

以下のfileは、選択routeに応じて設定directory以下に作られます。停止routeはmapper outputを作らず、BAM経路はcanonical FASTQ layerを必要としません。

metadataとdownload記録:

- **Raw ENA metadata**
  - Path: `<filereport-dir>/filereport_read_run_PRJNA<ID>_raw_tsv.txt`
  - ENA read_run metadataの生fileです。sample filter前にENAから何が返ってきたかを監査できます。
- **Selected ENA metadata**
  - Path: `<filereport-dir>/filereport_read_run_PRJNA<ID>_tsv.txt`
  - filter後のdownload scopeを定義するmetadataです。その後modality/platform routingでmapper scopeが絞られる場合は、対応assignment reportを確認してください。
- **Spreadsheet metadata**
  - Path: `<filereport-dir>/PRJNA<ID>.csv`
  - 選択metadataをspreadsheetで確認しやすくしたCSVです。
- **Platform inference JSON**
  - Path: `<filereport-dir>/platform_inference_PRJNA<ID>.json`
  - platform推定の記録です。metadata/FASTQそれぞれの判定、選択platform、evidence、halt判断を含みます。
- **GEO SOFT cache**
  - Path: `<filereport-dir>/geo_soft/`
  - metadata-based platform推定に使うGEO SOFT cacheです。
- **Generated download script**
  - Path: `<download-script-outputdir>/filereport_read_run_PRJNA<ID>_tsv_download_srr.sh`
  - 再現性確認や手動確認用に保存されるSRA download scriptです。
- **Project log**
  - Path: `<temporary-sra-download-dir>/prjna<ID>/prjna<ID>_log.txt`
  - project単位のlogです。download進行、推定message、warning、halt actionが残ります。

raw FASTQとread-structure記録:

- **Sample FASTQ directory**
  - Path: `<final-file-dir>/prjna<ID>/<sample_alias>/`
  - metadata filter後のsample別FASTQ directoryで、`SRRxxxx_1.fastq.gz` などのlocal名を使います。download済みfileはmapper preparation時にrename・書換えしませんが、archive basenameとadapterが付けたlocal名が同じとは限りません。
- **Read-structure assignment TSV**
  - Path: uniformなproject-level routeでは `<final-file-dir>/prjna<ID>/read_structure_assignment.tsv`、deferred non-10x routeでは `<final-file-dir>/prjna<ID>/<sample_alias>/read_structure_assignment.tsv`
  - source FASTQ suffixをcanonical logical role (`I1`, `I2`, `R1`, `R2`) に割り当てた記録で、scopeはrouteにより異なります。mapper preparationは該当fileを使い、raw FASTQをrenameせずcanonical symlinkを作ります。
- **Read-structure assignment JSON**
  - Path: uniformなproject-level routeでは `<final-file-dir>/prjna<ID>/read_structure_assignment.json`、deferred non-10x routeでは `<final-file-dir>/prjna<ID>/<sample_alias>/read_structure_inference.json`
  - 該当するassignmentとinference evidenceのJSON記録で、programmatic inspectionに使えます。
- **SRR FASTQ rearrangement manifest**
  - Path: `<final-file-dir>/prjna<ID>/srr_fastq_rearrangement.tsv`
  - source file名を保持したまま、downloadしたSRR FASTQをGSM/sample directoryへどのように割り当てたかを記録します。source/destinationのdevice、inode、size、mtime、ctime fingerprintも保存します。すでにFASTQ全体を検証済みで、同一filesystem内の移動によりdevice/inode、size、mtimeが保持された場合は、移動manifestとdestination ctimeを照合して、整理後にFASTQ全体を再展開せず以前の検証結果を安全に再利用します。filesystemをまたぐcopy、置換、fingerprint変化があればfull validationを再実行します。

BAM rescue input:

- **BAM rescue manifest**
  - Path: `<final-file-dir>/prjna<ID>/bam_inputs_manifest.tsv`
  - `--resolve-bam` で作られるmanifestです。downloadしたsubmitted BAM path、検出されたSAM tag、調査record数、raw barcode/UMI tagが完全だったrecord数を記録します。
- **Submitted BAM files**
  - Path: `<final-file-dir>/prjna<ID>/<sample_alias>/*.bam`
  - ENA/GEO source URLから直接downloadしたsubmitted BAM/alignment fileです。raw FASTQではありません。
- **BAM rescue summary**
  - Path: `<final-file-dir>/prjna<ID>/.uniscflow_bam_rescue.json`
  - project単位のBAM rescue summaryです。

mapper script、manifest、mapping output:

- **Mapper manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/mapper_inputs_manifest.tsv`
  - sampleごとのmanifestです。platform、実際に使うtarget mapper、要求されたtarget、生成script path、出力先を記録します。grouped Smart-seq2をSTARsolo SmartSeq manifestで安全に実装する場合など、上位の要求targetと実際のSTAR modeが異なる経路も監査できます。
- **Sample-modality assignment**
  - Path: `<mapper-output-dir>/prjna<ID>/sample_modality_assignment.tsv`
  - mixed-assayの各選択sampleについてmodality、GEX採用、明示的non-GEX除外、未解決の `manual_review` 判定とsample-level evidenceを記録します。scope判定の記録であり、mapper完了を示すものではありません。
- **Sample-platform routing manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/sample_platform_routing.tsv`
  - 解決済みmixed-platform projectについて、選択された各GSMのselected platform、endpoint、return code、reasonを1つずつ記録します。
- **Per-sample terminal endpoint**
  - Path: `<mapper-output-dir>/prjna<ID>/sample_route_endpoints/<GSM>.json`
  - sample-levelの `documented_halt`、`non_target_stop`、`unsupported_stop` routeで作られます。そのGSMにmapper commandを生成しなかった理由と、利用可能な場合はprofile固有guidanceを記録します。
- **STARsolo command**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/starsolo/command.sh`
  - 10x、Drop-seq、Seq-WellなどSTARsolo対応sampleの再現可能なSTARsolo commandです。
- **Canonical mapper FASTQs**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/<target>/fastqs/`
  - 生成されたmapper commandが読む軽量symlink layerです。STARsolo droplet/UMI targetではcanonical `R1` がbarcode/UMI、canonical `R2` がcDNAです。STAR + featureCountsやSalmon targetでは、canonical `R1/R2` がfull-length read pairまたはsingle-end readです。source FASTQはmapper preparationでcopyもoverwriteもされません。
- **Canonical FASTQ manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/<target>/fastqs/canonical_fastqs.tsv`
  - canonical FASTQ symlink、canonical role、lane番号、source role、source pathを記録します。
- **Manual-review source FASTQ manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/manual_review/fastqs/source/source_fastqs.tsv`
  - 低レベルのmapper-input helperをmanual-review profileで手動実行した時に、UniScFlowがsource FASTQを整理できるがmapper-readyな `R1/R2/I1/I2` roleを安全に割り当てられない場合に書かれることがあります。symlinkは `SRC1`、`SRC2` などのsource-role名を使い、元のdownload FASTQ名は `<final-file-dir>` 以下に保持されます。
- **Manual-review README**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/manual_review/README.txt`
  - 低レベルhelperで作成したmanual-review mapper-input directoryに付随することがあります。automatic mappingをhaltした理由と、後続のvendor-specific/manual preprocessing用にcanonicalまたはsource-role symlink FASTQを作成したかどうかを記録します。
- **10x sample-level inference record**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/starsolo/sample_level_10x_inference.json`
  - Cell Ranger chemistry定義とbarcode fileが渡された10x sampleで書かれます。sample-levelのchemistry call、barcode whitelist score、source role assignment、STARsolo command作成に使ったcanonical mapper FASTQ directoryを記録します。
- **STARsolo output**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/starsolo/starsolo_out/`
  - STAR log、`Solo.out`、raw/filtered matrix、barcode/UMI summaryが入ります。UniScFlowはSTARsoloを `--soloFeatures Gene GeneFull` で実行するため、STARが対応する場合は `Solo.out/Gene/` と `Solo.out/GeneFull/` の両方が作られます。
- **Mapper run manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/mapper_run_manifest.tsv`
  - mapper実行後に書かれます。sampleごとのcommand status、exit code、validation reasonを記録します。STARsolo、STAR + featureCounts、Salmon、generated Cell Ranger targetはそれぞれendpointを検証します。matrix targetでは空でないdimensionとfeature/barcode companion fileを要求し、exit 0だけでは成功扱いしません。
- **Mapping completion receipt**
  - Path: `<mapper-input-dir>/.uniscflow_mapping_complete.json`
  - 検証済みmapper outputをproject、selected SRR、filereport scope、platform、target、reference context、command checksumに結び付けます。UniScFlowはoutputを再検証し、このreceiptがcurrent run contextと一致する間だけ再利用します。
- **Smart-seq2 command**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/star_featurecounts/command.sh`
  - Smart-seq2/full-length sample用の再現可能なSTAR + featureCounts commandです。
- **Smart-seq2 output**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/star_featurecounts/star_featurecounts_out/`
  - Smart-seq2のSTAR alignmentとfeatureCounts output directoryです。
- **Smart-seq2 standardized matrix**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/star_featurecounts/star_featurecounts_out/uniscflow_matrix/`
  - UniScFlowで標準化したSmart-seq2 matrixです。`features.tsv`, `barcodes.tsv`, `matrix.mtx`, `counts.tsv` を含みます。STARsolo/10x風に、安定した統合keyとして `gene_id`、annotationとして `gene_name` を持ちます。
- **Grouped-sample manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/sample_groups.tsv`
  - `--sample-map-tsv` を指定した時に作られます。GSM/well/cell directoryをどの生物学的sampleにまとめたかを記録します。
- **Per-sample grouping manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_id>/sample_map.tsv`
  - 生物学的sampleごとの対応表です。
- **Grouped Smart-seq2 STARsolo manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_id>/mapper_inputs/starsolo/read_files_manifest.tsv`
  - well/cellごとのR1 FASTQ、R2 FASTQまたは `-`、cell IDを記録します。
- **Grouped Smart-seq2 command**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_id>/mapper_inputs/starsolo/command.sh`
  - 生物学的sampleごとに1回だけ実行する `STAR --soloType SmartSeq` commandです。
- **Grouped Smart-seq2 matrix**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_id>/mapper_inputs/starsolo/starsolo_out/Solo.out/Gene/` および `Solo.out/GeneFull/`
  - well/cell identityを保持したsample単位のgene x cell STARsolo matrixです。UniScFlowはSTARsoloを `--soloFeatures Gene GeneFull` で実行します。

STARsolo sampleとSTAR + featureCounts sampleでは、`--write-web-summary` により以下が作られます。

- **STARsolo web summary**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/starsolo/starsolo_out/web_summary.html`
  - STARsolo metric、mapping rate、barcode-rank plot、matrix statistics、実行commandをまとめた1ファイルHTML QC reportです。
- **STAR + featureCounts web summary**
  - Path: `<mapper-output-dir>/prjna<ID>/<sample_alias>/mapper_inputs/star_featurecounts/star_featurecounts_out/web_summary.html`
  - STAR alignment metric、featureCounts assignment category、detected genes、top expressed genes、count distribution、実行commandをまとめた1ファイルHTML QC reportです。
- **Web summary manifest**
  - Path: `<mapper-output-dir>/prjna<ID>/web_summary_manifest.tsv`
  - どのsampleでweb summaryを作ったか、どのmapper outputを使ったかを記録するmanifestです。
- **Multiplex metadata audit**
  - Path: `<mapper-output-dir>/prjna<ID>/multiplex_audit.json`
  - `confirmed`、`suspected`、`feature_companion`、`not_detected`、`unavailable`の判定と、その根拠になったmetadataを記録します。確定、疑い、feature companionの場合はrun logとSTARsolo GEX web summaryにもwarningを表示しますが、mapping endpointは変更しません。
- **Documented halt (recognized stop) web summary**
  - Path: `<mapper-output-dir>/prjna<ID>/halt_summary.html`
  - 現在のinput scopeに一致するhalt markerについて、profile固有の停止要因、必要な追加入力、推奨workflow、具体的な次の手順、再開条件を示すproject-level HTML reportです。

追加の監査file:

- **Inference report TSV**
  - Path: `--inference-report-tsv` で指定したpath
  - project-level 10x read-inference pathが作るoptional append-only TSVです。全platform共通のaudit logではなく、deferred non-10x routeは上記のsample別assignment fileを使います。
- **Download halt marker**
  - Path: `<final-file-dir>/prjna<ID>/.uniscflow_halt_after_download.json`
  - documented halt (recognized stop)または認識済みunsupported stopに対するscope-bound halt markerです。21のdocumented halt (recognized stop) profileでは、mapper command生成とmappingに進まなかった理由に加え、不足する入力、推奨workflow、具体的な次の手順、再開条件を記録します。unsupported assayでは `halt_type=unsupported_platform` または `halt_type=non_target_data` を記録し、profile resourceを与えればmapping可能になるとは主張しません。

## sampleの選び方

PRJNA内の全sampleを対象にする:

```bash
--ids 804521
```

1つのGSMだけ:

```bash
--sample-alias GSM5876625
```

複数GSM:

```bash
--sample-alias GSM5876624,GSM5876625
```

他のENA metadata列でfilter:

```bash
--filter run_accession=SRR17944814
--filter sample_title=Oc4mu
```

指定したGSMや値が見つからない場合は、以下のようなwarningを出します。

```text
Warning: No rows matched sample_alias=GSM0000000
```

## platform推定

`--platform auto` は、repository metadataをFASTQ構造・配列evidence、または適格なtag付きBAM evidenceと照合します。sample protocol、data processing、library field、title、characteristics、supplementary file名を調べ、判定根拠をinference reportに残します。

初期GEO summaryはdefaultで先頭3 GSMを使います（`--geo-soft-max-samples`）。sample-level routingでは、検証済みcacheと関連GSE family SOFTから選択GSM全体の証拠も収集し、不足sampleを記録します。default cacheは `<filereport-dir>/geo_soft` で、`--geo-soft-dir` で変更できます。cacheのintegrityを検証し、一時的取得失敗をretryし、GEOを利用できない場合はENA metadataをfallback evidenceとして使います。fallbackでも各routeの判定条件を満たす必要があります。

証拠の一致を基に自動routingし、見かけの矛盾はplatform固有の照合規則で評価します。例えば、十分な10x chemistry evidenceは曖昧なmetadataのchemistry名を解決できます。assay/platformの矛盾が解決しなければ、後から調べたsourceを無条件で優先せずreviewのため停止します。

Smart-seq2では、長いpaired-end FASTQだけからsingle-cell libraryとbulkまたはlow-input SMART-Seq RNA-seqを安全に区別できません。そのためUniScFlowは、automatic mappingへ進む前に追加のsingle-cell evidenceを要求します。ENA metadataの `TRANSCRIPTOMIC SINGLE CELL`、GSM/sample-levelのGEO本文中にあるsingle-cell RNA-seq、single-nucleus RNA-seq、FACS-sorted single cells、index sortingなどの明示的な記載をevidenceとして扱います。96/384-well plate processingの記載は、多数のGSM aliasまたはplate well座標を持つone-well-per-GSM型の構造を伴う場合にsingle-cell plate evidenceとして採用します。1つのGSEにsingle-cell libraryとbulk libraryが混在することがあるため、Series-levelのsingle-cell wordingだけでは自動mappingへ進みません。このようなsample-levelまたはplate-layout evidenceを伴わないSMART-Seq keywordだけでは、scRNA-seqとしてmappingせず意図的に停止します。

```bash
uniscflow --mode infer-platform \
  --ids 267857 \
  --platform auto \
  --filereport-dir /path/to/filereport \
  --final-file-dir /path/to/ready_fastqs
```

10xでは、metadataとFASTQの両方が10xを示していればplatformは一致とみなします。v2/v3などのchemistry/versionはmetadataよりもFASTQ whitelist/Cell Ranger chemistry推定を優先します。公共metadataではv2/v3表記がずれることがあるためです。

`--cellranger-chemistry-defs` と `--cellranger-barcodes-dir` がある場合、mapper preparationでもGSM/sampleごとに10x chemistryとread roleを再確認します。同一project内でbarcode/UMI streamのsuffixが異なる場合も、source FASTQをwhitelist・chemistry定義と照合し、canonical `R1` をbarcode/UMI、`R2` をcDNAとして準備します。download済みFASTQは上書きせず、対応を `sample_level_10x_inference.json` とcanonical manifestに記録します。

ChromiumとVisiumのprotocol記載が混在していても、独立したsample-local GEX証拠と全runのwhitelist・chemistry・read-role証拠が一致すれば、10x GEXとして解決できます。spatialの記載と解決時のwarningはinference reportに残し、whitelist一致だけで明示的なspatialやnon-GEX assayを上書きしません。[軽量10xデモ](docs/tutorials/quickstart_10x_prjna825585.md#read-the-inference-and-matrix-results)でこの例を紹介しています。

矛盾が解決しない場合は、sample-level protocol、選択run、read-role evidenceを先に確認してください。`--force-platform` はexpert向けoverrideであり、通常の矛盾解決手順ではありません。不足する実験resourceを補うことも、要求routeの生物学的妥当性を証明することもできません。

Drop-seqとSeq-Wellでは、`--platform dropseq` または `--platform seqwell` を指定すると、`--index1`, `--index2`, `--read1`, `--read2` を明示しない限りread構造推定が自動で有効になります。典型的な2-read layoutでは `_1` をbarcode/UMI `R1`、`_2` をcDNA `R2`、index readを `NULL` として扱います。明示的なplatform metadataがあり、選択した全runがこのcanonical layoutを満たす場合は、barcode readが20 bpより長くても有効です。CB12+UMI8 profileは先頭20塩基を使い、R1末尾の追加配列を無視します。UniScFlowは各runのread長とambiguous baseを独立にQCし、mapping thresholdは満たすもののwarning thresholdを下回るrunをreportに記録します。

FASTQ構造だけでも、「fixed whitelistを持たないdroplet UMI系」らしい、というfamilyまではかなり分かります。一方で、Drop-seqとSeq-Wellのように同じbarcode/UMI geometryを使うplatform名までは、FASTQだけでは安全に区別できません。そのためdefaultではmetadataで名前付きplatformを決めきれない場合は停止します。generic STARsolo routeを使う場合は、`--platform generic_droplet_umi` と完全なgeometryをユーザーが明示します。UniScFlowがFASTQ read長から座標を推測・最適化することはありません。すべての座標が正であること、cell barcodeとUMIが重ならず同じlogical read上にあること、cDNAがもう一方のlogical readにあること、各barcode FASTQでsampling readの70%以上がCB/UMI必要末端に達することを検証します。70%-<90%はwarningを残してmappingし、70%未満はhaltします。whitelist指定時のFASTQ別一致率warningはmapper profileとoptional HTML summaryへ残します。名前付きplatformを示す明示metadataまたは強いFASTQ evidenceと矛盾する場合はoverrideせず停止します。Drop-seqまたはSeq-Well metadataは、完全な明示geometryがその名前付きprofileと完全一致する場合だけcompatibleとして扱い、異なるgeometryならhaltします。DNBelab C4/C Series系は `dnbelab_c4` として認識し、DNBelab/PISA固有のpreprocessingが必要になり得るためdefaultではdownload + haltにします。

logical R1の先頭12 bpをcell barcode、続く8 bpをUMI、logical R2をcDNAとして明示する例:

```bash
uniscflow --mode all --ids PRJNA... \
  --platform generic_droplet_umi \
  --generic-cell-barcode-read R1 \
  --generic-cell-barcode-start 1 \
  --generic-cell-barcode-length 12 \
  --generic-umi-read R1 \
  --generic-umi-start 13 \
  --generic-umi-length 8 \
  --generic-cdna-read R2 \
  --star-index /path/to/star_index \
  --genes-gtf /path/to/genes.gtf
```

ここでの `R1`/`R2` はUniScFlowがcanonicalizeした後のlogical mapper roleです。barcode/UMI streamとcDNA streamがread長で一意ならsource suffixを推定できます。両方のsource streamが長いなど曖昧な場合は `--read1` と `--read2` も明示してください。generic geometryが曖昧なsource-file assignmentをoverrideすることはありません。

## Direct FASTQ fallback integrity

通常のSRA/ODP downloadが失敗した場合、UniScFlowはENA `fastq_ftp` やNCBI SDL source objectに登録された直接FASTQ.gzへfallbackできます。このdirect FASTQはdefaultで `--fastq-integrity-check gzip --fastq-integrity-retries 2` により検証されます。`gzip` checkは圧縮file全体を読むため、途中で切れたFASTQ.gzをread-structure inferenceへ流さず、削除して再downloadします。ENA `fastq_ftp` fileでは、ENAが値を提供している場合にdownload byte数とMD5を `fastq_bytes` / `fastq_md5` と照合します。

NCBI SDL fallbackには1回の実行あたり6時間のwall-clock上限があります（SDL lookupは1回90秒）。上限に達するとworkerを停止し、完了済み・途中までのfileを保持して、次回の実行では途中のfileから再開します。特に大きなsource fileでは、起動前に上限を延ばしてください。例:

```bash
export UNISCFLOW_SDL_STAGE_TIMEOUT_SECONDS=172800   # 48時間
```

`UNISCFLOW_SDL_LOOKUP_TIMEOUT_SECONDS` も同じ方法でlookup上限を調整できます。download stage全体は `docs/workflow.md` を参照してください。

FASTQ全体の検証結果はpath、size、mtime、ctimeとともにcacheします。sample別directoryへの整理では、内容が同一でもpathとctimeが変わります。このためUniScFlowは `srr_fastq_rearrangement.tsv` で移動前後のfingerprintを結び、同一inodeの安全な移動と確認できた場合だけ以前のfull validationを再利用します。coverage reportにはdirect cache hit、移動後cache hit、full validation件数を記録します。

## Submitted BAM Rescue

公共10x projectの一部では、raw split FASTQではなく `possorted_genome_bam.bam` や `*_alignments.bam` のようなalignment BAMがsubmitted fileとして登録されています。この場合、`fasterq-dump` から1本のFASTQ-like outputが得られても、元の10x R1/R2/I1構造は復元できません。

UniScFlowはdefaultでsubmitted BAMを確認し、利用可能なBAMをdownload・検証してSAM tagを検査します。適格なBAMにはraw barcode/UMIの配列・quality tagが必要です。有効なcurrent inputでcoverできないrunはSRA取得と欠落runのfallbackへ進みます。初期BAM検出を無効にするには `--no-resolve-bam` を指定します。

defaultでは `--bam-integrity-check full --bam-integrity-retries 2` が使われます。`full` はBAM全体を `samtools view -c` で読むため `samtools quickcheck` より遅いですが、STARsolo mapping中に初めて出るようなCRC/inflate errorも事前に検出できます。input coverage validationでも、BAMのpath、size、mtime、ctimeが変化した場合は同じcheckを再実行し、変化していないBAMだけintegrity cacheを再利用します。速度を優先する場合だけ `--bam-integrity-check quickcheck` を指定してください。

BAM rescue validationは、`samtools` が `PATH` 上で使える環境で実行してください。付属のconda環境とDocker imageには `samtools` が入っています。integrity validationはfail-closedで、custom環境に `samtools` がなければ `samtools_not_found` を記録し、未検証BAMを受理せず失敗させます。ENA metadataに値があるsubmitted BAMではbyte数とMD5も照合します。

Barcode/UMI tagはBAM全体で見つかったtag名の和集合ではなく、record単位で検証します。自動STARsolo BAM rescueには、調査した全recordで有効なraw `CR/UR`配列と同じ長さの文字列qualityが必要です。qualityは標準の`CY/UY`、または旧Cell Ranger 1.0-1.1の`CQ/UQ`に限ります。旧形式の受理には、STAR command内の旧referenceと`CELLRANGER_CS/CELLRANGER/EXTRACT_READS`入力由来の両方、および一様なbarcode/UMI長の記録が必要です。不正・部分的なcanonical qualityからaliasへはfallbackしません。これはmetadata-freeのplatform同定に必要な明示的`cellranger count` commandの代わりにはなりません。record単位のtag集計がないmanifestは自動rescueには使えないため、download/preparation stageを再実行してBAMを再検査し、`bam_inputs_manifest.tsv`を更新してください。

```text
CR/CY: raw cell barcode sequence/quality
UR/UY: raw UMI sequence/quality
CB/UB: corrected barcode/UMI。監査には有用ですが、自動raw BAM rescueには不十分です。
```

`CR/CY/UR/UY` が揃っている場合、生成されるSTARsolo commandは、headerを含む1本のSAM streamとして `--readFilesType SAM SE`, `--soloInputSAMattrBarcodeSeq CR UR`, `--soloInputSAMattrBarcodeQual CY UY` を指定します。submitted BAMが1本なら `samtools view -h -F 0x900` でstreamし、同一sampleに複数のrun-level BAMがある場合は `samtools merge` で1本に統合してから同じprimary-alignment filterを通します。これによりraw barcode/UMI tagを保ったままsecondaryおよびsupplementary alignmentを除外します。由来を確認した旧形式では代わりに`--soloInputSAMattrBarcodeQual CQ UQ`と記録済みbarcode/UMI長を指定し、選択BAM間でschemaと旧形式の長さが一致することを要求します。corrected `CB/UB`への置換やqualityの合成は行いません。このrouteの入力はraw FASTQではなく公開されたalignment BAMであり、input-run coverageとQCの判定は変更しません。

## STAR-based mapping

default workflowでは、まずSTAR genome indexを作ります。

```bash
uniscflow --mode build-star-index \
  --star-index /path/to/star-mm10-2020-A \
  --genome-fasta /path/to/refdata-gex-mm10-2020-A/fasta/genome.fa \
  --genes-gtf /path/to/refdata-gex-mm10-2020-A/genes/genes.gtf \
  --threads 30
```

GTFはSTAR index作成時に使います。STARsolo mapping時には `--sjdbGTFfile` を再度渡さず、`--genomeDir` 内の `geneInfo.tab`, `transcriptInfo.tab` などをSTARsoloに読ませます。Smart-seq2のSTAR + featureCounts mappingでは、featureCounts実行時に `--genes-gtf` も使います。

次に、platformごとのSTARベースmappingを実行します。10x、Drop-seq、Seq-WellではSTARsolo commandを生成します。10x datasetの例:

```bash
uniscflow --mode mapping \
  --ids 804521 \
  --platform 10x \
  --filereport-dir /path/to/filereport \
  --download-script-outputdir /path/to/download_script \
  --temporary-sra-download-dir /path/to/sra_tmp \
  --final-file-dir /path/to/ready_fastqs \
  --mapper-output-dir work/mapper_ready \
  --star-index /path/to/star-mm10-2020-A \
  --genes-gtf /path/to/refdata-gex-mm10-2020-A/genes/genes.gtf \
  --cellranger-chemistry-defs /path/to/cellranger/chemistry_defs.json \
  --cellranger-barcodes-dir /path/to/cellranger/barcodes \
  --threads 30
```

10x datasetでは、Cell Ranger chemistry定義とbarcode directoryが使える場合、`--starsolo-whitelist` を通常は指定しなくて大丈夫です。UniScFlowが推定chemistryに対応するwhitelistをSTARsolo用に自動選択します。手動で固定したい場合だけ `--starsolo-whitelist` を指定します。

`--mode mapping` はsample別mapper directoryを準備し、`command.sh` を実行します。`--dry-run` はworkflow stageのcommandを表示します。実際のmapper commandを実行前に確認するには、`--mode prepare` で生成した `command.sh` を読んでください。

Drop-seqやSeq-Wellの例:

```bash
uniscflow --mode mapping \
  --ids PRJNAxxxxxx \
  --platform dropseq \
  --filereport-dir /path/to/filereport \
  --download-script-outputdir /path/to/download_script \
  --temporary-sra-download-dir /path/to/sra_tmp \
  --final-file-dir /path/to/ready_fastqs \
  --mapper-output-dir work/mapper_ready \
  --star-index /path/to/star-mm10-2020-A \
  --genes-gtf /path/to/refdata-gex-mm10-2020-A/genes/genes.gtf \
  --threads 16
```

複数sampleを同時に走らせたい場合は、CPUとmemoryに余裕がある時だけ `--run-mapper-parallel N` を指定します。

## STAR-based web summary

mapping後に、Cell Ranger風の1ファイルHTML QC reportをsampleごとに生成できます。STARsolo outputではbarcode/UMI中心のsummaryを作り、Smart-seq2のSTAR + featureCounts outputでは `Log.final.out`, `counts.txt`, `counts.txt.summary` からfull-length RNA-seq用のsummaryを作ります。

```bash
uniscflow --mode report \
  --ids PRJNAxxxxxx \
  --mapper-output-dir work/mapper_ready
```

出力先は各sampleのmapper output directoryです。

```text
work/mapper_ready/prjnaxxxxxx/GSMxxxx/mapper_inputs/starsolo/starsolo_out/web_summary.html
work/mapper_ready/prjnaxxxxxx/GSMxxxx/mapper_inputs/star_featurecounts/star_featurecounts_out/web_summary.html
work/mapper_ready/prjnaxxxxxx/web_summary_manifest.tsv
```

STARsolo reportは `Log.final.out`, `Solo.out/Gene/Summary.csv`, `UMIperCellSorted.txt`, `Gene` filtered matrix を読みます。`GeneFull` は下流解析用に生成されますが、web summaryのdefault matrixには使いません。STAR + featureCounts reportは `Log.final.out`, `featurecounts/counts.txt`, `featurecounts/counts.txt.summary` を読みます。sample-level metrics、plot、output statistics、実行したcommandをまとめます。

![STARsolo web summary preview](docs/assets/starsolo-web-summary-preview.png)

例: [PRJNA1087433 / GSM8146611 のSTARsolo web summary](docs/examples/starsolo_web_summary_PRJNA1087433_GSM8146611.html)。

mapping後に自動で作る場合は `--write-web-summary` を付けます。

```bash
uniscflow --mode mapping \
  --ids PRJNAxxxxxx \
  --platform dropseq \
  --filereport-dir /path/to/filereport \
  --download-script-outputdir /path/to/download_script \
  --temporary-sra-download-dir /path/to/sra_tmp \
  --final-file-dir /path/to/ready_fastqs \
  --mapper-output-dir work/mapper_ready \
  --star-index /path/to/star-mm10-2020-A \
  --genes-gtf /path/to/refdata-gex-mm10-2020-A/genes/genes.gtf \
  --threads 16 \
  --write-web-summary
```

## Mapper script生成

input取得・推定後、`--mode prepare` でmappingせずmapper inputとscriptを生成できます。FASTQ経路はcanonical link、適格なtag付きBAM経路はSAM streamを使います。current halt markerがあるscopeでは通常の準備は停止します。resource fileを追加するだけで解除されると考えず、記録された次の手順に従ってください。

```bash
uniscflow --mode prepare \
  --ids 804521 \
  --platform 10x \
  --target starsolo \
  --filereport-dir /path/to/filereport \
  --download-script-outputdir /path/to/download_script \
  --temporary-sra-download-dir /path/to/sra_tmp \
  --final-file-dir /path/to/mapper_ready_fastqs \
  --mapper-output-dir work/mapper_ready \
  --star-index /path/to/star_index \
  --genes-gtf /path/to/genes.gtf \
  --cellranger-chemistry-defs /path/to/cellranger/chemistry_defs.json \
  --cellranger-barcodes-dir /path/to/cellranger/barcodes \
  --threads 30
```

Drop-seqやSeq-WellではSTARsolo用のcommand templateを作ります。

```bash
uniscflow --mode prepare \
  --ids PRJNAxxxxxx \
  --platform dropseq \
  --target starsolo \
  --filereport-dir /path/to/filereport \
  --download-script-outputdir /path/to/download_script \
  --temporary-sra-download-dir /path/to/sra_tmp \
  --final-file-dir /path/to/ready_fastqs \
  --mapper-output-dir work/mapper_ready \
  --star-index /path/to/star_index \
  --genes-gtf /path/to/genes.gtf \
  --threads 16
```

標準CLIでは、実行可能なmapping scriptを生成する前に、実在する `--star-index` と、それを作成した時に使用したものと対応する `--genes-gtf` の両方が必要です。`uniscflow --mode check` はfull run前に両方のpathを検証します。

Smart-seq2ではdefaultのSTAR + featureCounts経路を要求します。granularity検証により、GSM-as-cellはSTAR + featureCounts、解決済みrun-as-cellやreview済みcell groupはSTARsolo SmartSeqへroutingされます。したがって以下の例は、counting backendを無条件に固定するものではありません。

```bash
uniscflow --mode prepare \
  --ids PRJNAxxxxxx \
  --platform smartseq2 \
  --target star_featurecounts \
  --filereport-dir /path/to/filereport \
  --download-script-outputdir /path/to/download_script \
  --temporary-sra-download-dir /path/to/sra_tmp \
  --final-file-dir /path/to/ready_fastqs \
  --mapper-output-dir work/mapper_ready \
  --star-index /path/to/star_index \
  --genes-gtf /path/to/genes.gtf \
  --threads 16
```

STAR + featureCountsが選ばれた場合、scriptはGSM/cellごとに `STAR` を実行し、`Aligned.sortedByCoord.out.bam` を作ってから `featureCounts` を実行します。paired-end inputは `featureCounts -p --countReadPairs` によりread pair/fragment単位でcountし、single-endはpaired-end optionを使いません。

featureCounts終了後、UniScFlowは `star_featurecounts_out/uniscflow_matrix/` を作り、Smart-seq2 outputもSTARsolo/10x風のgene row conventionにそろえます。

```text
features.tsv  # gene_id, gene_name, Gene Expression
barcodes.tsv  # GSM-as-cell経路のGSM/cell識別子
matrix.mtx    # gene x cell count matrix（各GSM outputは1列）
counts.tsv    # gene_id/gene_name/count のtable表示
```

platformをまたいだ統合では、`gene_id` をprimary keyにし、`gene_name` は表示用annotationとして扱います。gene symbolは変わったり重複したりすることがある一方、同一reference内のgene IDは安定したrow identifierとして使いやすいためです。

### 1つのGSMが1 wellを表すplate project

Smart-seq2/full-length projectで1 GSMが1 well/cellと判定できる場合、sample mapなしでSTAR + featureCountsによる細胞単位のmappingを行います。複数GSMを指定した生物学的sample単位の行列にまとめたい場合は、任意で `--sample-map-tsv` を利用できます。UniScFlowはsampleごとにSTARsolo SmartSeq manifestとmapping commandを1つずつ作り、各GSM/well/cellを別々の行列列として保持します。

これはSTAR公式のplate-based Smart-seq mode [`--soloType SmartSeq --readFilesManifest`](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#plate-based-smart-seq-scrna-seq) を使う実装です。

明確にrun-as-cellと判定できるdepositもuser mapなしでSTARsolo SmartSeqへroutingできます。細胞単位の判定と生物学的sampleのgroupingは別であり、生物学的な所属は自動推定しません。cell granularityが曖昧な場合は推測せず停止します。

TSVでは、GSM列に `gsm_accession`, `gsm`, `sample_alias`, `geo_accession`、生物学的sample列に `sample_id`, `biological_sample`, `sample_name`, `group_id`, `sample` のいずれかを使います。任意のcell/well列は `cell_id`, `well_id`, `well`, `cell`, `library_id`、condition列は `condition`, `treatment`, `group` です。

`sample_map.tsv` の例:

```tsv
gsm_accession	sample_id	cell_id	condition
GSM000001	treated_mouse_1	A01	treated
GSM000002	treated_mouse_1	A02	treated
GSM000003	control_mouse_1	A01	control
```

自分の確認済みmapを使う場合は、以下のplaceholder accessionを置き換え、選択GSMを網羅する `sample_map.tsv` を用意してください。具体的な3細胞の例は、[Smart-seq2 tutorialの任意grouping節](docs/tutorials/quickstart_smartseq2_prjna701252.md#optional-group-cells-with-a-sample-map)を参照してください。

```bash
uniscflow --mode all \
  --ids PRJNAxxxxxx \
  --platform auto \
  --sample-map-tsv sample_map.tsv \
  --filereport-dir work/filereport \
  --download-script-outputdir work/download_script \
  --temporary-sra-download-dir work/tmp_sra \
  --final-file-dir work/raw \
  --mapper-output-dir work/mapper_ready \
  --star-index /path/to/star_index \
  --genes-gtf /path/to/genes.gtf \
  --threads 16 \
  --write-web-summary
```

上の形式例にあるgroup名を使う場合、出力は以下のようにまとまります。

```text
work/mapper_ready/prjna<ID>/
  sample_groups.tsv
  treated_mouse_1/
    sample_map.tsv
    mapper_inputs/starsolo/
      read_files_manifest.tsv
      command.sh
      starsolo_out/Solo.out/Gene/
      starsolo_out/Solo.out/GeneFull/
  control_mouse_1/
    sample_map.tsv
    mapper_inputs/starsolo/
      read_files_manifest.tsv
      command.sh
      starsolo_out/Solo.out/Gene/
      starsolo_out/Solo.out/GeneFull/
```

`--sample-map-tsv` を指定した場合、選択されたすべてのGSM directoryが対応表に含まれている必要があります。これにより、生物学的sample単位のGSMと、1 well/cell単位のGSMを誤って混ぜることを防ぎます。

## read構造の指定と推定

手動で指定する場合:

```bash
--index1 1 \
--index2 NULL \
--read1 2 \
--read2 3
```

これは以下を意味します。

```text
*_1.fastq.gz -> I1
*_2.fastq.gz -> R1
*_3.fastq.gz -> R2
```

10xでchemistry-awareに自動推定する場合:

```bash
--cellranger-chemistry-defs /path/to/cellranger/lib/python/cellranger/chemistry_defs.json \
--cellranger-barcodes-dir /path/to/cellranger/lib/python/cellranger/barcodes \
--min-barcode-match-rate 0.8
```

chemistryが分かっている場合は候補を絞れます。

```bash
--cellranger-chemistry SC3Pv3-polyA
```

推定結果をTSVに残す:

```bash
--inference-report-tsv work/read_inference.tsv
```

## Optional legacy Cell Ranger mapping

Cell Rangerが入っている既存のDocker containerを指定してmappingします。このlegacy経路は10x Genomicsと推定されたprojectだけに限定されます。current filereportのsample/SRRだけを選び、記録済みread assignmentから一時的な相対canonical FASTQ symlinkを作り、実行後に除去します。filtered gene-expression matrixがない場合は成功扱いしません。失敗sampleは `cellranger_run_manifest.tsv` に記録し、後続sampleは続行しますが、project全体はdefaultでnon-zero終了します。部分outputを明示的に許容する場合だけ `--allow-partial-success` を使ってください。`--cellranger-include-introns true|false` でCell Ranger versionごとのdefaultを上書きできます。

### 動作確認済みのCell Ranger 7.1.0 container

動作確認済みの構成では、Docker Hubの `nfcore/cellranger:7.1.0` imageを使います。

```bash
docker pull nfcore/cellranger:7.1.0

docker run -itd \
  --name cellranger710 \
  --hostname cellranger710 \
  -v /path/to/cellranger710_workspace:/home/yard \
  nfcore/cellranger:7.1.0 \
  /bin/bash
```

このmountでは、host側の:

```text
/path/to/cellranger710_workspace
```

がCell Ranger container内では:

```text
/home/yard
```

として見えます。

例:

```text
/path/to/cellranger710_workspace/raw  ->  /home/yard/raw
/path/to/cellranger710_workspace/ref  ->  /home/yard/ref
```

Cell RangerおよびDocker imageは、それぞれのlicense/distribution termsに従って使ってください。

```bash
uniscflow --mode mapping \
  --ids 804521 \
  --platform auto \
  --filereport-dir /path/to/filereport \
  --download-script-outputdir /path/to/download_script \
  --temporary-sra-download-dir /path/to/sra_tmp \
  --final-file-dir /path/to/cellranger710_workspace/raw \
  --mapping-engine cellranger \
  --cellranger-container cellranger710 \
  --transcriptome /home/yard/ref/refdata-gex-mm10-2020-A \
  --localcores 30 \
  --localmem 500 \
  --dir-in-container /home/yard/raw \
  --file-dir-in-host /path/to/cellranger710_workspace/raw \
  --dir-in-host /path/to/cellranger710_workspace/raw \
  --download-source SRR \
  --no-bam
```

mapping前にmount対応を確認します。

```bash
uniscflow --mode check \
  --ids 804521 \
  --mapping-engine cellranger \
  --filereport-dir /path/to/filereport \
  --download-script-outputdir /path/to/download_script \
  --temporary-sra-download-dir /path/to/sra_tmp \
  --final-file-dir /path/to/cellranger710_workspace/raw \
  --cellranger-container cellranger710 \
  --transcriptome /home/yard/ref/refdata-gex-mm10-2020-A \
  --dir-in-container /home/yard/raw \
  --file-dir-in-host /path/to/cellranger710_workspace/raw \
  --dir-in-host /path/to/cellranger710_workspace/raw
```

標準UniScFlow imageにはDocker CLIが入っていないため、Docker socketをmountするだけでは動きません。legacy Cell Ranger mappingはhostから実行するか、互換Docker clientを追加した派生imageをbuildし、その場合にsocketをmountしてください。

```bash
-v /var/run/docker.sock:/var/run/docker.sock
```

## validation

```bash
uniscflow --mode validate \
  --ids 804521 \
  --filereport-dir work/filereport \
  --download-script-outputdir work/download_script \
  --final-file-dir work/raw
```

`--mode validate` はinputの準備状態を確認します。

- sampleが0件になっていないか
- download scriptが空ではないか
- raw-input directoryが存在し、metadataとoutput pathが対応するか
- 選択runのintegrity確認済みcoverage。該当する場合は検証済みreceiptやscopeに対応したhalt evidenceも扱います

このmodeはmappingを実行せず、mapper runnerが行うtarget固有のmatrix検証の代わりにはなりません。

## config file

基本的にはcommand-line optionで指定できます。繰り返し使う環境ではconfig fileも使えます。同梱exampleの値はすでにdefaultなので、install方法にかかわらず以下を使えます。

```bash
uniscflow --mode plan
```

repository checkout内ではexampleを明示指定できます。checkout外では、自分で用意したTOMLのabsolute pathを指定してください。

```bash
uniscflow --mode plan --config ./config/example.toml
```

command-line optionはconfig fileの値を上書きします。

## repository構成

```text
bin/uniscflow                           command-line entrypoint
uniscflow.py                            workflow controller
tools/legacy/                           現行runtime helper（directory名は歴史的なもの）
docs/workflow.md                        workflow詳細
docs/tutorials/                         user向けtutorial
docs/examples/                          実行可能なtutorial script
tests/                                  softwareの自動test
docker/cellranger-v10-check/             Cell Ranger chemistry file確認用container
Dockerfile                              workflow Docker image
environment.yml                         conda environment
```

## 注意点

- Cell Ranger本体は含めません。
- barcode whitelistは再配布しません。Cell Ranger installationなど、利用可能なsourceから取得してください。
- STAR/STARsoloとfeatureCountsはconda環境とDocker image内でBiocondaからinstallします。詳細は `THIRD_PARTY_NOTICES.md` を参照してください。
- manifest依存platformでは、外部barcode、well、sample manifestが必要です。UniScFlowはFASTQだけから無理にmappingせず、意図的に停止します。
- 公共archiveのmetadataは不整合があるため、大規模download前に選択sampleを確認してください。
- 自動read推定は、新しいchemistryや特殊なdatasetでは必ず確認してください。
- uniform remapping自体はstudy間統合、batch/生物学的効果の除去、donor/replicateの同定を行いません。後続解析前にcell qualityとgroupingを確認してください。
- 大規模実行ではsource download、中間file、canonical linkの参照元、mapper outputを置く容量が必要です。linkを利用中はsource FASTQを保持してください。

## Citation

[凍結済みUniScFlow 1.0.0 release](https://doi.org/10.5281/zenodo.22271248)とUniScFlow論文を引用してください。citation metadataは [`CITATION.cff`](CITATION.cff)、release手順は [`docs/release.md`](docs/release.md) にあります。このDOIは保存された版を識別し、以後の `main` の変更は含みません。解析で使った正確なGit commitとenvironmentまたはcontainer digestを記録してください。

## License

UniScFlow独自のsource codeは[BSD 3-Clause License](LICENSE)で提供されます。third-party tool、reference resource、barcode file、container componentにはそれぞれのlicenseが適用されます。詳細は[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)を参照してください。
