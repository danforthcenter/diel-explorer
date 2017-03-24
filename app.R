library(shiny)
library(ggplot2)
library(grid)
library(shinydashboard)
library(reshape)
library(reshape2)
library(stringr)
library(heatmaply)
library(gridExtra)
library(gplots)
library(svglite)

########################################################################################
# Read Setaria Circadian Data In (See methods for details on generation of datatables)
# See setaria-circadian.R script for how data was cleaned up / normalized
########################################################################################

ldhhf.llhcf<-read.table(file='data/setaria-ldhhf.llhcf.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)
ldhhf.llhcf.norm<-read.table(file='data/setaria-ldhhf.llhcf.norm.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)
entrainment<-read.table(file='data/entrainment.information.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)

########################################################################################
# Setaria Shiny Application -ui.R
########################################################################################
ui <- dashboardPage(
  
  dashboardHeader(title="Diel Explorer"),

  dashboardSidebar(disable = T),

  dashboardBody(
    
    tabsetPanel(
      
      ########################################################################################
      tabPanel(title="Welcome",
           fluidRow(
             box(title="Diel Explorer",width=12, solidHeader = T,status = 'primary',
                 p("This tool is brought to you by the",a(href="http://www.gehanlab.org/",target='_blank',"Gehan Lab"),
                   "at the Donald Danforth Plant Science Center. For the code used to generate this app,
                   please visit ", a(href="https://github.com/maliagehan/diel-explorer/",target='_blank',"our Github"),". For more information on how the data was processed please refer
                   to XXXXX(Pub/Preprint). To download the raw data go here:"),
                 
                 actionButton("rawdata","Go to Raw Data",icon=icon("th"),
                              onclick ="window.open('http://google.com', '_blank')")
                 
                ),
             
             box(title="Using this Tool",width=12, solidHeader = T,status = 'primary',
                  p("In the SAMPLE INFO tab see the available datasets and conditions"),
                  hr(),
                  p("In the SEARCH AND BROWSE DATA tab above, you can search for genes of interest using  by either
                  the search bar or by uploading a .txt file with gene ids. Alternatively, you can use the data filters
                  to browse the data.Once you have loaded or searched for your selections the data can be viewed and dowloaded. 
                  A plot of the data can be seen on the Plot Data tab once the plot button is hit.
                  Warning: asking to plot a large line graph is messy, we would suggest switching to a heatmap"),
                  hr(),
                  p("In the PLOT DATA tab above, data selected from Select Data or Browse Data tabs is plotted"),
                  hr(),
                  p("In the ADDING YOUR OWN DATA tab above, we briefly describe how to alter this Shiny App so you
                    can graph your own timecourse data.")
                )
               )
              ),
      ########################################################################################
      tabPanel(title="Sample Info",
               fluidRow(
                 box(title="Sample Info",width=12, solidHeader = T,status = 'danger',
                     p("Notes for Reference XXX: Setaria viridis was grown for 10 days then released into constant conditions for circadian sampling
                     every 2 hours for 48 hours. Again, for more details sampling conditions, please refer to XXXXX(Pub/Preprint).
                     The two entrainment conditions for this experimental set are as follows:"),
                     dataTableOutput("conditions.table")
                 )
               )
      ),
      ########################################################################################
      tabPanel(title="Search and Browse Data",
               fluidRow(
                 box(title="Search Data with GENEID or GO",width=12, solidHeader = T,status = 'success',collapsible = TRUE, collapsed = TRUE,
                     h4("Search using small sets of GENEIDs"),
                     p("GENEIDs, Orthologs, or GO separated by a comma are allowed"),
                     textInput('gene_search_text', label="example: Sevir.2G310200.1, Sevir.1G000100.1", value = ""),
                     h4("Search using small sets of GO TERMS"),
                     textInput('go_search_text', label="example:GO:GO:0008270", value = ""),
                     h4("Search using small sets of ORTHOLOG GENEIDs"),
                     textInput('orth_search_text', label="example:AT3G17930.1,LOC_Os01g59080.1", value = ""),
                     p("refresh page to clear search")
                 )),
               
               fluidRow(
                 box(title="Search Data with File",width=12, solidHeader = T,status = 'success',collapsible = TRUE,collapsed=TRUE,
                     h4("Upload a file of GENEIDs, the GENEIDs should not be quoted, one line per geneid"),
                     fileInput('file.geneid', 'Choose file to upload',
                               accept = c('text/csv','text/comma-separated-values','text/tab-separated-values','text/plain', '.csv','.tsv')),
                     checkboxInput('header', 'Header', FALSE),
                     h4("Upload a file of ORTHOLOG GENEIDs, the ORTHOLOG GENEIDs should not be quoted, one line per geneid"),
                     fileInput('file.ortholog', 'Choose file to upload',
                               accept = c('text/csv','text/comma-separated-values','text/tab-separated-values','text/plain', '.csv','.tsv')),
                     checkboxInput('header1', 'Header', FALSE),
                     p("refresh page to clear search")
                 )),
               fluidRow(
                 box(title="Genes, Orthologs, or GO Selected with Search",width=12, solidHeader = T,status = 'success',collapsible = TRUE,collapsed =TRUE,
                     verbatimTextOutput("selected.genes"),
                     verbatimTextOutput("selected.orthologs"),
                     verbatimTextOutput("selected.go")
                     )),
               fluidRow(
                 box(title="Browse and Filter Data",width=12, solidHeader = T,status = 'success',
                     fluidRow(
                            column(3,
                                   radioButtons('data.norm', 'Normalize Data',c(Yes='TRUE',NO='FALSE'),'TRUE')),
                            column(3,
                                    selectInput("species1",
                                                "Species:",
                                                c("All", 
                                                 sort(unique(as.character(ldhhf.llhcf$species)))))),
                            column(3,
                                   selectInput("dataset1",
                                               "Entrainment:",
                                               c("All", 
                                                 sort(unique(as.character(ldhhf.llhcf$dataset)))))),
                            column(3,
                                   selectInput("bhq.cutoff1",
                                               "Benjamini-Hochberg Q-Value:",
                                               c("All",sort(unique((ldhhf.llhcf$bhq.cutoff)))))),
                            column(3,
                                   selectInput("adjp.cutoff1",
                                               "Adjusted P-Value:",
                                               c("All",sort(unique((ldhhf.llhcf$adjp.cutoff)))))),
                            column(3,
                                   selectInput("period1",
                                               "Period:",
                                               c("All",
                                                 sort(unique(as.character(ldhhf.llhcf$PERIOD)))))),
                            column(3,
                                   selectInput("lag1",
                                               "Lag (Phase):",
                                               c("All",
                                                 sort(unique(as.character(ldhhf.llhcf$LAG)))))),
                            
                            column(3,
                                   downloadButton('download.selected', "Download Selected Data"))
                            ),
                          div(style = 'overflow-x: scroll', dataTableOutput("setaria.data"))
                  )
                 )
                ),
      ########################################################################################
      tabPanel(title="Plot Data",
               fluidRow(
                 box(title="Plot Data",width=12, solidHeader = T,status = 'info',
                     h3("Warning: Attempting to plot too much data on line graph will be slow/unresponsive and messy."),
                     actionButton("plot.data",label="Plot Selected Data as Line Graph"),
                     actionButton("plot.heat", label="Plot Selected Data as Heatmap" ),
                     radioButtons('colorby', 'Color Line Graph By',c(Dataset='TRUE',GENEID='FALSE'),'TRUE'),
                     radioButtons('rowcol', 'Scale Heatmap',c(Row='TRUE',Column='FALSE'),'TRUE'),
                     radioButtons('average', 'Average Replicates',c(Yes='TRUE',No='FALSE'),'TRUE'),
                     downloadButton('download.plot', "Download Line Graph"),
                     downloadButton('download.heat',"Download Heatmap"),
                     plotOutput("circadian.expression.plot"),
                     plotlyOutput("circadian.expression.heat")
                )
               )
             ),
      ########################################################################################
      tabPanel(title="Adding Your Own Data",
               fluidRow(
                 box(title="Adding Your Own Data",width=12, solidHeader = T,status = 'warning',collapsible = TRUE,
                     p(" We encourage you to use this frame work to explore your own data and hope that you will add it to the
                       public repository as well. If you have any problems adding your own data or with the exisiting data,
                       please ", a(href="https://github.com/maliagehan/diel-explorer/issues",target='_blank',"contact us"),". The following are 
                       the basic steps to add data.")
                 )),
               fluidRow(
                 box(title="STEP 1:Data Format",width=12, solidHeader = T,status = 'warning',collapsible = TRUE,
                     p("This takes in the resulting output of ",a(href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3119870/",target='_blank',"JTK Cycle"),
                       ". For more information on the parameters used to process data with JTK Cycle, please refer to XXX (preprint/pub).
                       The resulting .txt files are then further processed in Step 2.")
                 )),
               fluidRow(
                 box(title="STEP2: Add Data and Format Further",width=12, solidHeader = T,status = 'warning',collapsible = TRUE,
                     p("Your sampling frequencing and headers may differ from ours. So you can alter your naming scheme to match ours, 
                       (if it makes sense to) or you might need to replace your data with ours. Replacing our data with yours is a bit more involved, because,
                       it will change the header names for the graphing steps downstream. For now we will assume that your data headers are like ours.
                       There is a R script in the data folder called setaria-circadian.R. The setaria-circadian.R script was used to add some additional information to the JTK Cycle output. 
                       For example, two JTK Cycle outputs were joined together with an rbind step then merged with an annotation file. 
                       We add in columns for species name, categories for BH.Q and ADJ.P, and then find the maximum expression
                       value per replicate so that normalized data can be calculated. Once all of these things are done, we have have the raw (expression values not normalized to maximum)
                       and relative/normalized expression data tables that are used as input data. Once you add the necessary columns your dataset, you could rbind it to to the Setaria dataset
                       and reload the app.")
                 )),
               fluidRow(
                 box(title="STEP3: Add Sample Info",width=12, solidHeader = T,status = 'warning',collapsible = TRUE,
                     p("There is a file called entrainment.information.txt. Add a row to that file with your sample information and save it.")
                 ))
            ),
      ########################################################################################
      tabPanel(title="Contact Us",
               fluidRow(
                 box(title="Contact US",width=12, solidHeader = T,status = 'danger',
                     p("For questions or if you are interested in adding your own data contact ",
                       a(href="https://github.com/maliagehan/diel-explorer/issues",target='_blank',"Malia Gehan")),
                     p("For mor information on the Gehan lab visit our ",
                       a(href="http://www.gehan-lab.org/",target='_blank',"website"))
                 )
               )
            )
    )
  )
)
    
########################################################################################
# Setaria Shiny Application -server.R
########################################################################################
server<-function(input,output){
  
  #output for sample info tab #########################################################
  
  output$conditions.table<-renderDataTable(entrainment, options=list(paging=FALSE,searching=FALSE))
  
  #output for select data tab #########################################################
  
  #output for browse data tab #########################################################
  
  #get the search terms
  searchgenes<-reactive({
    genes<-gsub(" ","",input$gene_search_text)
    genes1<-data.frame(strsplit(genes,","))
    colnames(genes1)<-c("GENEID")
    genes1
    })
  
  #get the file contents
  searchfile<-reactive({
    geneids<-input$file.geneid
    
    if(is.null(geneids))
      return(NULL)
    
    genes1<-read.csv(geneids$datapath,header=input$header, strip.white = TRUE)
    colnames(genes1)<-c("GENEID")
    genes1
  })
  
  #join the gene search and the file contents
  joinedsearch<-reactive({
    rbind(searchgenes(),searchfile())
  })
  
  #output the list of genes so the user can see it
  output$selected.genes<-renderPrint({
    joinedsearch()
  })
  
  #get the search ortholog terms
  searchorthologs<-reactive({
    orth<-gsub(" ","",input$orth_search_text)
    orth1<-data.frame(strsplit(orth,","))
    colnames(orth1)<-c("ORTHOLOGS")
    orth1
  })
  
  #get the ortholog file terms
  searchfileorth<-reactive({
    orth<-input$file.ortholog
    
    if(is.null(orth))
      return(NULL)
    
    orth1<-read.csv(orth$datapath,header=input$header1,strip.white =TRUE)
    colnames(orth1)<-c("ORTHOLOGS")
    orth1
  })
  
  #join the ortholog search and the ortholog file contents
  joinedsearchorth<-reactive({
    rbind(searchorthologs(),searchfileorth())
  })
  
  #output the list of genes so the user can see it
  output$selected.orthologs<-renderPrint({
    joinedsearchorth()
  })
  
  #get list of go terms
  searchgo<-reactive({
    go<-gsub(" " ,"",input$go_search_text)
    go1<-data.frame(strsplit(go,","))
    colnames(go1)<-c("GO-TERM")
    go1
  })
  
  #output the list of go so the user can see it
  output$selected.go<-renderPrint({
    searchgo()
  })
  
  #filter data
  setaria.input<-reactive({
    
    setaria.table<-ldhhf.llhcf.norm

    if(input$data.norm=='FALSE'){
      setaria.table<-ldhhf.llhcf
    }
    
    if(nrow(joinedsearch())!=0 | nrow(searchgo())!=0 | nrow(joinedsearchorth())!=0){
       colnum<-ncol(setaria.table)
       setaria.subset <- data.frame(matrix(ncol=colnum, nrow = 0))
       colnames(setaria.subset) <- paste0(c(colnames(setaria.table)))
    
       if(nrow(joinedsearch())!=0){
         for(x in 1:nrow(joinedsearch())){
           search<-as.character(joinedsearch()[x,1])
           row<-setaria.table[(grep(search,setaria.table$GENEID)),]
           setaria.subset<-rbind(row,setaria.subset)
         }
       }
       
       if(nrow(joinedsearchorth())!=0){
         for(x in 1:nrow(joinedsearchorth())){
           search<-as.character(joinedsearchorth()[x,1])
           row<-setaria.table[(grep(search,setaria.table$ortholog)),]
           setaria.subset<-rbind(row,setaria.subset)
         }
       }
    
       if(nrow(searchgo())!=0){
         for(x in 1:nrow(searchgo())){
           search<-as.character(searchgo()[x,1])
           row<-setaria.table[(grep(search,setaria.table$GO)),]
           setaria.subset<-rbind(row,setaria.subset)
         }
       }
       
       setaria.table<-setaria.subset
    }

    if(input$species1 !="All"){
      setaria.table<-setaria.table[setaria.table$species==input$species1,]
    }
    
    if(input$dataset1 !="All"){
      setaria.table<-setaria.table[setaria.table$dataset==input$dataset1,]
    }
    
    
    if(input$bhq.cutoff1!="All"){
      setaria.table<-setaria.table[setaria.table$BH.Q<=as.numeric(input$bhq.cutoff1),]
    }
    
    if(input$adjp.cutoff1!="All"){
      setaria.table<-setaria.table[setaria.table$ADJ.P<=as.numeric(input$adjp.cutoff1),]
    }
    
    if(input$period1 !="All"){
      setaria.table<-setaria.table[setaria.table$PERIOD==input$period1,]
    }
    
    if(input$lag1 !="All"){
      setaria.table<-setaria.table[setaria.table$LAG==input$lag1,]
    }
    setaria.table
  })
    
  output$setaria.data<-renderDataTable(({
    setaria.input()}),options=list(searching=FALSE))
  
  output$download.selected <- downloadHandler(
    filename = function() { paste('circadian_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),'.csv', sep='') },
    content = function(file) {
      write.csv(setaria.input(), file)
      }
    )

  #output for plot data tab ###########################################################  
  
  plot.input1<-reactive({
    
           setaria.plot1<-setaria.input()
           setaria.plot1$uniqueid<-NA
           setaria.plot1$uniqueid<-paste(setaria.plot1$GENEID,"-",setaria.plot1$dataset,sep="")
           
           setaria.rep1<-subset(setaria.plot1,select=c(GENEID,uniqueid,dataset,ZT02_rep1,ZT04_rep1,ZT06_rep1,ZT08_rep1,ZT10_rep1,ZT12_rep1,
                                                             ZT14_rep1,ZT16_rep1,ZT18_rep1,ZT20_rep1,ZT22_rep1,ZT24_rep1,ZT26_rep1,
                                                             ZT28_rep1,ZT30_rep1,ZT32_rep1,ZT34_rep1,ZT36_rep1,ZT38_rep1,ZT40_rep1,
                                                           ZT42_rep1,ZT44_rep1,ZT46_rep1,ZT48_rep1))
           
           colnames(setaria.rep1)<-c("GENEID","uniqueid","dataset","ZT02","ZT04","ZT06","ZT08","ZT10","ZT12",
                                     "ZT14","ZT16","ZT18","ZT20","ZT22","ZT24","ZT26",
                                     "ZT28","ZT30","ZT32","ZT34","ZT36","ZT38","ZT40",
                                     "ZT42","ZT44","ZT46","ZT48")
           setaria.rep1$rep<-"rep1"
           
           setaria.rep2<-subset(setaria.plot1,select=c(GENEID,uniqueid,dataset,ZT02_rep2,ZT04_rep2,ZT06_rep2,ZT08_rep2,ZT10_rep2,ZT12_rep2,
                                                       ZT14_rep2,ZT16_rep2,ZT18_rep2,ZT20_rep2,ZT22_rep2,ZT24_rep2,ZT26_rep2,
                                                       ZT28_rep2,ZT30_rep2,ZT32_rep2,ZT34_rep2,ZT36_rep2,ZT38_rep2,ZT40_rep2,
                                                       ZT42_rep2,ZT44_rep2,ZT46_rep2,ZT48_rep2))
           
           colnames(setaria.rep2)<-c("GENEID","uniqueid","dataset","ZT02","ZT04","ZT06","ZT08","ZT10","ZT12",
                                     "ZT14","ZT16","ZT18","ZT20","ZT22","ZT24","ZT26",
                                     "ZT28","ZT30","ZT32","ZT34","ZT36","ZT38","ZT40",
                                     "ZT42","ZT44","ZT46","ZT48")
           setaria.rep2$rep<-"rep2"
           
           setaria.rep1.melt<-melt(setaria.rep1,by=c(uniqueid))
           setaria.rep2.melt<-melt(setaria.rep2,by=c(uniqueid))
           
           if(input$average=='TRUE'){
             setaria.all<-rbind(setaria.rep1,setaria.rep2)
             setaria.all.melt<-melt(setaria.all,by=c(uniqueid,rep))
             setaria.all.average<-melt(cast(setaria.all.melt,uniqueid~variable,mean))
             
             setaria.all.average$dataset<-NA
             setaria.all.average$GENEID<-NA
             
             gid<-str_split_fixed(setaria.all.average$uniqueid,"-",2)
             
             setaria.all.average$dataset<-gid[,2]
             setaria.all.average$GENEID<-gid[,1]
             
             if(input$colorby=='FALSE'){
               ggplot()+
                 geom_point(data=setaria.rep1.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(uniqueid), shape="REP1"))+
                 geom_point(data=setaria.rep2.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(uniqueid), shape="REP2"))+
                 geom_line(data=setaria.all.average,aes(x=variable,y=value,group=factor(uniqueid),color=factor(uniqueid)))+
                 theme_bw()+
                 labs(title="Circadian Expression", y="EXPRESSION or RELATIVE EXPRESSION", x="TIME(ZT)")+
                 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right")}else{
               ggplot()+
                 geom_point(data=setaria.rep1.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(dataset), shape="REP1"))+
                 geom_point(data=setaria.rep2.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(dataset), shape="REP2"))+
                 geom_point(data=setaria.all.average,aes(x=variable,y=value,group=factor(uniqueid),color=factor(dataset),shape=factor(GENEID)))+
                 geom_line(data=setaria.all.average,aes(x=variable,y=value,group=factor(uniqueid),color=factor(dataset)))+
                 theme_bw()+
                 labs(title="Circadian Expression", y="EXPRESSION or RELATIVE EXPRESSION", x="TIME(ZT)")+
                 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right")}
             }else{
           if(input$colorby=='FALSE'){
             ggplot()+
               geom_point(data=setaria.rep1.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(GENEID), shape="REP1"))+
               geom_point(data=setaria.rep2.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(GENEID), shape="REP2"))+
               geom_line(data=setaria.rep1.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(GENEID)))+
               geom_line(data=setaria.rep2.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(GENEID)))+
               theme_bw()+
               labs(title="Circadian Expression", y="EXPRESSION or RELATIVE EXPRESSION", x="TIME(ZT)")+
               theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right")
           }else{
             ggplot()+
               geom_point(data=setaria.rep1.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(dataset), shape="REP1"))+
               geom_point(data=setaria.rep2.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(dataset), shape="REP2"))+
               geom_line(data=setaria.rep1.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(dataset)))+
               geom_line(data=setaria.rep2.melt,aes(x=variable,y=value,group=factor(uniqueid),color=factor(dataset)))+
               theme_bw()+
               labs(title="Circadian Expression", y="EXPRESSION or RELATIVE EXPRESSION", x="TIME(ZT)")+
               theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right")
           }}
           })
  
  plot.input2<-reactive({
    
    setaria.plot1<-setaria.input()
    setaria.plot1$uniqueid<-NA
    setaria.plot1$uniqueid<-paste(setaria.plot1$GENEID,"-",setaria.plot1$dataset,sep="")
    
    setaria.rep1<-subset(setaria.plot1,select=c(GENEID,uniqueid,dataset,ZT02_rep1,ZT04_rep1,ZT06_rep1,ZT08_rep1,ZT10_rep1,ZT12_rep1,
                                                ZT14_rep1,ZT16_rep1,ZT18_rep1,ZT20_rep1,ZT22_rep1,ZT24_rep1,ZT26_rep1,
                                                ZT28_rep1,ZT30_rep1,ZT32_rep1,ZT34_rep1,ZT36_rep1,ZT38_rep1,ZT40_rep1,
                                                ZT42_rep1,ZT44_rep1,ZT46_rep1,ZT48_rep1))
    
    colnames(setaria.rep1)<-c("GENEID","uniqueid","dataset","ZT02","ZT04","ZT06","ZT08","ZT10","ZT12",
                              "ZT14","ZT16","ZT18","ZT20","ZT22","ZT24","ZT26",
                              "ZT28","ZT30","ZT32","ZT34","ZT36","ZT38","ZT40",
                              "ZT42","ZT44","ZT46","ZT48")
    setaria.rep1$rep<-"rep1"
    
    setaria.rep2<-subset(setaria.plot1,select=c(GENEID,uniqueid,dataset,ZT02_rep2,ZT04_rep2,ZT06_rep2,ZT08_rep2,ZT10_rep2,ZT12_rep2,
                                                ZT14_rep2,ZT16_rep2,ZT18_rep2,ZT20_rep2,ZT22_rep2,ZT24_rep2,ZT26_rep2,
                                                ZT28_rep2,ZT30_rep2,ZT32_rep2,ZT34_rep2,ZT36_rep2,ZT38_rep2,ZT40_rep2,
                                                ZT42_rep2,ZT44_rep2,ZT46_rep2,ZT48_rep2))
    
    colnames(setaria.rep2)<-c("GENEID","uniqueid","dataset","ZT02","ZT04","ZT06","ZT08","ZT10","ZT12",
                              "ZT14","ZT16","ZT18","ZT20","ZT22","ZT24","ZT26",
                              "ZT28","ZT30","ZT32","ZT34","ZT36","ZT38","ZT40",
                              "ZT42","ZT44","ZT46","ZT48")
    setaria.rep2$rep<-"rep2"
    
    setaria.rep1.melt<-melt(setaria.rep1,by=c(uniqueid))
    setaria.rep2.melt<-melt(setaria.rep2,by=c(uniqueid))
    
    if(input$average=='TRUE'){
      setaria.all<-rbind(setaria.rep1,setaria.rep2)
      setaria.all.melt<-melt(setaria.all,by=c(uniqueid,rep))
      setaria.all.average<-melt(cast(setaria.all.melt,uniqueid~variable,mean))
      setaria.all.average<-dcast(data = setaria.all.average,formula = uniqueid~variable,fun.aggregate = sum,value.var = "value")
      setaria.all.average1<-subset(setaria.all.average,select=-c(uniqueid))
      row.names(setaria.all.average1)<-setaria.all.average$uniqueid
      setaria.all.matrix<-setaria.all.average1
      
      if(input$rowcol=='FALSE'){
          color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
          heatmaply(setaria.all.matrix, scale="column",dendrogram = "none",margins=c(40, 200),grid_color = "grey", colors=color.palette(1001))}
      else{
        color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
        heatmaply(setaria.all.matrix, scale="row",dendrogram = "none",margins=c(40, 200),grid_color = "grey",colors=color.palette(1001))
    }}
    else{
      setaria.all<-rbind(setaria.rep1,setaria.rep2)
      setaria.all$uniqueid<-paste(setaria.all$uniqueid,"-",setaria.all$rep,sep="")
      setaria.all<-setaria.all[order(setaria.all$uniqueid),]
      setaria.all.matrix<-subset(setaria.all,select=-c(GENEID,uniqueid,dataset,rep))
      row.names(setaria.all.matrix)<-setaria.all$uniqueid
      
      if(input$rowcol=='FALSE'){
        color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
        heatmaply(setaria.all.matrix, scale="column",dendrogram = "none",margins=c(40, 200),grid_color = "grey", colors=color.palette(256))
        }else{
        color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
        heatmaply(setaria.all.matrix, scale="row",dendrogram = "none",margins=c(40, 200),grid_color = "grey",colors=color.palette(256))
      }} 

  })
  
  plot.input3<-reactive({
    
    setaria.plot1<-setaria.input()
    setaria.plot1$uniqueid<-NA
    setaria.plot1$uniqueid<-paste(setaria.plot1$GENEID,"-",setaria.plot1$dataset,sep="")
    
    setaria.rep1<-subset(setaria.plot1,select=c(GENEID,uniqueid,dataset,ZT02_rep1,ZT04_rep1,ZT06_rep1,ZT08_rep1,ZT10_rep1,ZT12_rep1,
                                                ZT14_rep1,ZT16_rep1,ZT18_rep1,ZT20_rep1,ZT22_rep1,ZT24_rep1,ZT26_rep1,
                                                ZT28_rep1,ZT30_rep1,ZT32_rep1,ZT34_rep1,ZT36_rep1,ZT38_rep1,ZT40_rep1,
                                                ZT42_rep1,ZT44_rep1,ZT46_rep1,ZT48_rep1))
    
    colnames(setaria.rep1)<-c("GENEID","uniqueid","dataset","ZT02","ZT04","ZT06","ZT08","ZT10","ZT12",
                              "ZT14","ZT16","ZT18","ZT20","ZT22","ZT24","ZT26",
                              "ZT28","ZT30","ZT32","ZT34","ZT36","ZT38","ZT40",
                              "ZT42","ZT44","ZT46","ZT48")
    setaria.rep1$rep<-"rep1"
    
    setaria.rep2<-subset(setaria.plot1,select=c(GENEID,uniqueid,dataset,ZT02_rep2,ZT04_rep2,ZT06_rep2,ZT08_rep2,ZT10_rep2,ZT12_rep2,
                                                ZT14_rep2,ZT16_rep2,ZT18_rep2,ZT20_rep2,ZT22_rep2,ZT24_rep2,ZT26_rep2,
                                                ZT28_rep2,ZT30_rep2,ZT32_rep2,ZT34_rep2,ZT36_rep2,ZT38_rep2,ZT40_rep2,
                                                ZT42_rep2,ZT44_rep2,ZT46_rep2,ZT48_rep2))
    
    colnames(setaria.rep2)<-c("GENEID","uniqueid","dataset","ZT02","ZT04","ZT06","ZT08","ZT10","ZT12",
                              "ZT14","ZT16","ZT18","ZT20","ZT22","ZT24","ZT26",
                              "ZT28","ZT30","ZT32","ZT34","ZT36","ZT38","ZT40",
                              "ZT42","ZT44","ZT46","ZT48")
    setaria.rep2$rep<-"rep2"
    
    setaria.rep1.melt<-melt(setaria.rep1,by=c(uniqueid))
    setaria.rep2.melt<-melt(setaria.rep2,by=c(uniqueid))
    
    if(input$average=='TRUE'){
      setaria.all<-rbind(setaria.rep1,setaria.rep2)
      setaria.all.melt<-melt(setaria.all,by=c(uniqueid,rep))
      setaria.all.average<-melt(cast(setaria.all.melt,uniqueid~variable,mean))
      setaria.all.average<-dcast(data = setaria.all.average,formula = uniqueid~variable,fun.aggregate = sum,value.var = "value")
      setaria.all.average1<-as.matrix(subset(setaria.all.average,select=-c(uniqueid)))
      rownames(setaria.all.average1)<-setaria.all.average$uniqueid
      setaria.all.matrix<-setaria.all.average1
      return(setaria.all.matrix)}
    else{
      setaria.all<-rbind(setaria.rep1,setaria.rep2)
      setaria.all$uniqueid<-paste(setaria.all$uniqueid,"-",setaria.all$rep,sep="")
      setaria.all<-setaria.all[order(setaria.all$uniqueid),]
      setaria.all.matrix<-as.matrix(subset(setaria.all,select=-c(GENEID,uniqueid,dataset,rep)))
      rownames(setaria.all.matrix)<-setaria.all$uniqueid
      return(setaria.all.matrix)}
    })

  observeEvent(input$plot.data,{
    output$circadian.expression.plot<-renderPlot({
      grid.arrange(plot.input1(),ncol=1)
    })})
  
  observeEvent(input$plot.heat,{
    output$circadian.expression.heat<-renderPlotly({
      plot.input2()
    })})
  
  output$download.plot <- downloadHandler(
    filename = function() { paste('circadian_plot_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),".svg", sep='') },
    content=function(file){
      ggsave(file, plot = plot.input1(), device = "svg",height =8, width=10)
    }, contentType='image/svg')
  
  output$download.heat <- downloadHandler(
    filename = function() { paste('circadian_plot_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),".pdf", sep='') },
    content=function(file){
      color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
      if(input$rowcol=='FALSE'){scale1="column"}else{scale1="row"}
      pdf(file, width=8, height=7, useDingbats =FALSE)
      heatmap.2(plot.input3(),
                Rowv=FALSE,
                Colv=FALSE,
                dendrogram='none',
                scale=scale1,
                col=color.palette(256),
                trace='none',
                margins=c(8,20),
                symbreaks=FALSE,
                symm=FALSE,
                cex.main=0.75,
                density.info="none"
                )
      dev.off()
    },contentType='image/pdf')
  
  #output for add data tab ############################################################  

}

########################################################################################
#  Setaria Shiny Application - run app.
########################################################################################

# Run the application 
shinyApp(ui = ui, server = server)

########################################################################################
