
options(browser="firefox")
#.libPaths("/usr/lib/R/library")

#wor_dire="C:\\Users\\prunotto\\Desktop\\voynich"
#setwd(wor_dire)

library(edgeR)
library(combinat)
library(reshape2)
library(pheatmap)
library(matrixTests)
library(netstat)
library(RSelenium)
library(stringr)
library(stringi)
library(splitstackshape)
library(reshape2)
library(factoextra)
library(pheatmap)
library(ComplexHeatmap)
library(ComplexUpset)
library(dplyr)
library(tidytext)
library(bestNormalize)
library(ggrepel)
library(ggfortify)
library(readr)
library(showtext)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ggthemes)
library(Cairo)
library(ggpubr)

set.seed(1234)
font_add("voynich_font", "eva1.ttf")#"voynich-1.23-webfont.ttf")
showtext_auto()


counts_df=read.csv("voynich_counts.csv")

counts_ma=acast(ngrams~I,value.var="n",data=counts_df)
index=is.na(counts_ma)
counts_ma[index]=0

counts_df$S=dplyr::case_when(
	counts_df$I=="H"~"Herbal",
	counts_df$I=="P"~"Pharma",
	counts_df$I=="S"~"Stars",
	counts_df$I=="C"~"Cosmos",
	counts_df$I=="Z"~"Zodiac",
	counts_df$I=="B"~"Biology",
	counts_df$I=="T"~"Text",
	counts_df$I=="A"~"Astro",
	)



ups_plot=function(counts_df,var_sele,min_freq,min_class_freq,min_sect,max_show_word,ws){

	var_all=c("Astro","Cosmos","Stars","Zodiac","Biology","Herbal","Pharma","Text")

	index0=counts_df$n>=min_freq&counts_df$S%in%var_sele
	in_sections=aggregate(S~ngrams,data=counts_df,length)
	index1=counts_df$ngrams%in%in_sections[in_sections$S>=min_sect,]$ngrams
	index=index0&index1
	tdf=counts_df[index,]
	tdf_ma=acast(ngrams~S,value.var="n",data=tdf)
	index=is.na(tdf_ma)
	tdf_ma[index]=0
	tdf_ma=tdf_ma[rowSums(tdf_ma)>0,colSums(tdf_ma)>0]
	matr=as.data.frame(tdf_ma>0)

	matr=matr[,var_all[!is.na(var_all[match(var_all,colnames(matr))])]]


	base_upse=ComplexUpset::upset(
		matr,sort_sets=FALSE,
		sort_intersections_by=c('degree', 'cardinality'),
		colnames(matr)[colnames(matr)%in%var_sele],
		min_size=min_class_freq,
		min_degree=min_sect)


	temp=base_upse$patches$plots[[2]]$data
	bw=cbind.data.frame(word=rownames(temp),inte=temp$exclusive_intersection,ctrl=temp$exclusive_intersection_size)
	bw=bw[!is.na(bw$ctrl),]


	abw=list()
	lhm=list()
	for(i in unique(bw$inte))
	{
		bwi=bw[bw$inte==i,]$word
		bwc=tdf[tdf$ngrams%in%bwi,]
		th=aggregate(n~ngrams+S,data=bwc,sum)
		lh=acast(ngrams~S,value.var="n",data=th)
		lhm[[i]]=lh
		ah=aggregate(n~ngrams,data=bwc,sum)
		ah=tail(ah[order(ah$n),],max_show_word)
		abw[[i]]=ah
		abw[[i]]$upset_label=abw[[i]]$ngrams
	}
	wsel=unique(do.call(rbind,abw))
	af=aggregate(n~ngrams,data=tdf[tdf$ngrams%in%bw$word,c("ngrams","n")],sum)
	mf=merge(af,wsel,by=c("n","ngrams"),all.x=TRUE)

	toupse=merge(matr,mf,by.x="row.names",by.y="ngrams",all.x=TRUE)
	index=is.na(toupse$upset_label)
	toupse$upset_label[index]=""


	tit=paste("Words appearing in VM at least ",min_freq," times, in at least ",min_sect," subjects.",sep="")
	tim=paste("Combinations among at least ",min_sect," subjects, each containing at least ",min_class_freq, " words.",sep="")

	df=aggregate(word~inte,data=bw,length)
	df=df[match(rle(as.character(base_upse$data$intersection))$values,df$inte),]
	df$x=nrow(df):1

	upse=NULL
	upse=ComplexUpset::upset(
		toupse,
		colnames(toupse)[colnames(toupse)%in%var_sele],
		height_ratio=1,
		sort_sets=FALSE,
		sort_intersections_by=c('degree', 'cardinality'),
		wrap=TRUE,
		min_size=min_class_freq,
		min_degree=min_sect,
		set_size=FALSE,
		matrix=intersection_matrix()+ggtitle(tim),
		themes=upset_modify_themes(
			list(
				'intersections_matrix'=theme(
					axis.text.y=element_text(size=18),
					axis.title.x.bottom=element_text(size=0)
					)#,
				#'Word occurrence in VM'=theme(
				#	axis.text.y=element_text(size=15),
				#	#axis.title.y.right=element_text(size=25),
				#	axis.text.x=element_text(size=0),
				#	axis.title.x.bottom=element_text(size=0)					
				#	#axis.title.y=element_text(size=0)
				#	)
				)
			),
		base_annotations=list(),
		annotations=list(
			'Word occurrence in VM'=(
				ggplot(
					mapping=aes(y=n))+
				geom_jitter(
					position = position_jitter(seed = 1),
					aes(color=n), na.rm=TRUE,show.legend=FALSE,alpha=0.3)+
				geom_boxplot(alpha=0.3, na.rm=TRUE,show.legend=FALSE,outlier.shape=NA,col="grey")+
				geom_text_repel(
					position = position_jitter(seed = 1),
					aes(label=upset_label,family="voynich_font"),
					col="black",
					size=ws,alpha=0.9,max.overlaps=100,show.legend=FALSE)+	
				geom_text(data=df,inherit.aes=FALSE,y=-Inf,aes(x=x,label=word),vjust="inward",size=4)+
				scale_y_log10()+labs(title=tit)+
				rotate_y_text(angle=90)+
				theme(axis.text.y=element_text(size=12),axis.title.y.left=element_text(size=16))#+theme_bw()
				)
			)
    )#+ggtitle('Combinations among VM subjects')

	return(upse)

}




ui <- dashboardPage(
  dashboardHeader(title = "Voynich UpSet"),
  dashboardSidebar(),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    fluidRow(
    	column(2,
    		    pickerInput(
    				inputId = "S",
    				label ="Select subjects:", 
    				choices = unique(counts_df$S),
    				selected = unique(counts_df$S),
    				multiple=TRUE,
    				options = pickerOptions(actionsBox = TRUE)
    			),
    			sliderInput(
    				inputId= "min_freq",
    				label="Words occur at least # times:",
    				min=1,max=100,
    				value=4,
    				step=1
    				),
    			sliderInput(
    				inputId= "min_sect",
    				label="Words occur in at least # subjects:",
    				min=1,max=8,
    				value=4,
    				step=1
    				),
    			sliderInput(
    				inputId= "min_class_freq",
    				label="Minimal # words WITHIN UpSet class:",
    				min=1,max=100,
    				value=5,
    				step=1
    				), 
    			sliderInput(
    				inputId= "max_show_word",
    				label="Max # words SHOWN in UpSet class:",
    				min=0,
    				max=100,
    				value=1,
    				step=1
    				),     				   			
    			sliderInput(
    				inputId= "word_size",
    				label="Word size:",
    				min=1,max=10,
    				value=5,
    				step=0.5
    				)   			
    	),                            	
      column(10,
      	plotOutput(
      		"plot2", 
			dblclick = "plot2_dblclick",
    		brush = brushOpts(
    			id = "plot2_brush",
    			resetOnNew = TRUE),
      		height = 460)      	
      	)
    )
  )
)

server <- function(input, output) {

	  output$plot2=renderPlot({
	    var_sele=input$S
	    min_freq=input$min_freq
	    min_class_freq=input$min_class_freq	    
	    min_sect=input$min_sect
		max_show_word=input$max_show_word
	    ws=input$word_size
	    ups_plot(counts_df,var_sele,min_freq,min_class_freq,min_sect,max_show_word,ws)
	  })


}

voynich_app=shinyApp(ui, server)
voynich_app

