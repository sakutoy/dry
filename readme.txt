RNA-Seqの納品生データ（*fastq.gz）からDEG検出を行うパイプラインをまとめたシェルスクリプト

RUN環境の確認
OS：mac OS Catalina v10.15.7
シェル：zsh    *bashだとrunしないかも、要確認
requirments: python3, R, miniconda, trim-galore, hisat2, stringtie, samtools, deeptools
詳細はhttps://qiita.com/sakutoy/items/8c6461d675244aba0a94およびhttps://qiita.com/sakutoy/items/f6dcf8ff75f5dea53b56を参照

sample:今回はメダカの肝臓のbulk_RNA-Seqデータを使用する、cab1とcab5という2つのサンプルをペアエンドでシーケンスしている。

初期状態の確認
まずmacのデスクトップにanalysisというフォルダを作り、その中に以下のファイルがあるか確認する
・サンプルごとのペアエンドシーケンスデータ（今回はcab_liver_1_1.fastq.gzとcab_liver_1_2.fastq.gzおよびcab_liver_5_1.fastq.gzとcab_liver_5_2.fastq.gz）
・ゲノムリファレンスのファイル（Oryzias_latipes.ASM223467v1.dna.toplevel.fa）
・アノテーション情報を含むgtfファイル（Oryzias_latipes.ASM223467v1.104.gtf）
・prepDE.pyというファイル、なければインストールする（analysisへ移動してcurl -O https://ccb.jhu.edu/software/stringtie/dl/prepDE.pyを実行）
・パイプライン実行ファイル（DEG_analysis.sh）
・RでDEG検出＆プロット実行ファイル（deg_detect.R）、このファイルはDEG_analysis.shが呼び出して使う
・このファイル

バックアップ
・解析でエラーが出ると元データが失われてしまう可能性がある。必ず解析を始める前にanalysisのディレクトリを丸ごとどこか別の場所へコピペしてバックアップしておく

解析を始める
・初期状態の確認とバックアップが済んだら、ターミナルを開いてDEG_analysis.shを実行する。permission denied というエラーが出たら、chmod 755 [DEG_analysis.shのパス]を実行して権限を付与
・解析は数時間ほどかかる。うまくRUNすれば、degというフォルダにcsvファイル（gene_count_matrix.csvとtranscript_count_matrix.csv）が生成され、これらを用いてedgeRでDEG検出を行った結果（hoge.txtとプロット図）が出力される。

