#!/bin/sh

dir=`pwd`
cd "${dir}/Desktop/analysis"

sample1="cab_liver_1"
sample2="cab_liver_5"
sample_size=2
sample_array=($sample1 $sample2)

#sample1に個体1を、sample2に個体5のラベルを格納、配列にもサンプル情報を入れておく


#1.1 トリミングとクオリティチェック
for i in ${sample_array[@]};do
  echo $i
  #for文で配列から順に個体ラベル名を取り出して
  left=${i}_1.fastq.gz
  #left,rightにペアエンドのファイル名を格納
  right=${i}_2.fastq.gz
  trim_galore --q 20 --illumina --fastqc --trim1 --paired ${left} ${right}
  #それぞれのサンプルに対してペアエンドのリードをトリム＆QCを行う
  mkdir $i
  find . -type f -name "*${i}*" | xargs -I% mv % $i/.
  #cab1の生データと出力結果をcab1フォルダに格納する
done

#1.2　hisat2用リファレンスインデックスの作成
hisat2-build Oryzias_latipes.ASM223467v1.dna.toplevel.fa ASMol_index
mkdir index
find . -name '*ASMol*' | xargs -I% mv % index/.

#1.3 マッピング
for sample in ${sample_array[@]};do       #for文で配列から順に個体ラベル名を取り出して
  hisat2 -p 8 -x index/ASMol_index -1 ${sample}/${sample}_1_val_1.fq.gz -2 ${sample}/${sample}_2_val_2.fq.gz -S ${sample}_out.sam 2> ${sample}_out_res.txt
  #マッピングを実行
  samtools sort -@ 8 -O bam -o ${sample}_out.sort.bam ${sample}_out.sam
  mkdir ${sample}_hisat2_out
  find . -type f -name "*${sample}_out*" | xargs -I% mv % ${sample}_hisat2_out/.
done


#1.4 発現量取得
mkdir "ballgown"
for sample in ${sample_array[@]}
#for文で配列から順に個体ラベル名を取り出して
do
  stringtie ${sample}_hisat2_out/${sample}_out.sort.bam -e -B -p 4 -G Oryzias_latipes.ASM223467v1.104.gtf -o ballgown/${sample}/${sample}.gtf
done

#1.5 ファイル変換
mv prepDE.py ballgown
cd ballgown 　　　　　　　　　　　　　　　　　　　　　　　　
#curl -O https://ccb.jhu.edu/software/stringtie/dl/prepDE.py       #もしダウンロードできていなかったらコメントオンにする
2to3 -w prepDE.py
python3 prepDE.py

mkdir deg
#degというディレクトリを作成、deg解析はここで行う

cp ./ballgown/gene_count_matrix.csv ./deg/gene_count_matrix.csv
cp ./ballgown/transcript_count_matrix.csv ./deg/transcript_count_matrix.csv
cp ./deg_detect.R ./deg/deg_detect.R
#csvファイルとRスクリプトをdegディレクトリにコピーする

cd deg

Rscript deg_detect.R
