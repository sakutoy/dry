#このファイルのRUNの初期状態はhisat2およびstringtieでマッピングリードをカウントして発現量を取得し、gene_count_matrix.csvというファイルに記録した状態
#ディレクトリ構造はデスクトップ直下のanalysisにこのスクリプトとgene_count_matrix.csvが置かれている状態（DEG_analysisでそのように操作されるはず）




in_f <- "gene_count_matrix.csv"       #入力ファイル名を指定してin_fに格納
out_f <- "DEG_result.csv"                   #出力ファイル名を指定してout_fに格納
param_G1 <- 1                          #G1群のサンプル数を指定
param_G2 <- 1                          #G2群のサンプル数を指定
param_FDR <- 0.05                      #false discovery rate (FDR)閾値を指定

#必要なパッケージをロード
library(edgeR)                         #パッケージの読み込み

#入力ファイルの読み込みとラベル情報の作成
data <- read.table(in_f, header=TRUE, row.names=1, sep=",", quote="")#in_fで指定したファイルの読み込み
data <- as.matrix(data)                #データの型をmatrixにしている
data.cl <- c(rep(1, param_G1), rep(2, param_G2))#G1群を1、G2群を2としたベクトルdata.clを作成

#本番

d <- DGEList(counts=data,group=data.cl)#DGEListオブジェクトを作成してdに格納
d <- calcNormFactors(d)                #TMM正規化を実行
d <- estimateGLMCommonDisp(d, method="deviance", robust=TRUE, subset=NULL, verbose=TRUE)    #モデル構築(common Dispersionを算出)

out <- exactTest(d)                    #exact test (正確確率検定)で発現変動遺伝子を計算した結果をoutに格納
#tmp <- topTags(out, n=nrow(data), sort.by="none")#検定結果を抽出
p.value <- out$table$PValue            #p値をp.valueに格納
q.value <- p.adjust(p.value, method="BH")#q値をq.valueに格納
ranking <- rank(p.value)               #p.valueでランキングした結果をrankingに格納
sum(q.value < param_FDR)               #FDR閾値(q.value < param_FDR)を満たす遺伝子数を表示

#ファイルに保存(テキストファイル)
tmp <- cbind(rownames(data), data, p.value, q.value, ranking)#入力データの右側にp.value、q.value、rankingを結合した結果をtmpに格納
write.table(tmp, out_f, sep=",", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存

#AUC値を計算(高いほどよい方法; ただの検証用)  このへんはよくわからん
library(ROC)                           #パッケージの読み込み
param_DEG <- 1:sum(q.value < param_FDR)                    #DEGの位置を指定
obj <- rep(0, nrow(data))              #初期値として全てが0の(non-DEGに相当)ベクトルobjを作成
obj[param_DEG] <- 1                    #DEGの位置に1を代入
AUC(rocdemo.sca(truth=obj, data=-ranking))#AUC計算


#次にDEG（q.valueがFDR=0.05以下の転写産物）を識別して MAplotを描いてみる （参照：https://kazumaxneo.hatenablog.com/entry/2017/07/11/114021）
table<- as.data.frame(tmp)  　　　　　　　　　　　#データフレーム型に変換
is.DEG <- as.logical(table$q.value<param_FDR)　　　　　　#tableからq.valueの行を取り出して0.05より低いものをis.DEGに収納
DEG.names<-rownames(table)[is.DEG]　　　　　　　　#is.DEGの名前を取得
png("MA-plot.png")
plotSmear(out, de.tags=DEG.names, cex=0.3)　　　　　#DEGとそれ以外を識別してMAplotを出力
dev.off()
