# load libraries , scripts ,data

source("simulation.R")

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(GenomicRanges)
library(grid)
library(gridExtra)


rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r = rf(16)

mc = detectCores(all.tests = FALSE, logical = FALSE)

shinyServer(function(input,output){

  # make user specific variables
  reads <- reactive({

    # creates the GRanges object with the 5' positions    
    nreads = input$nreads
    S = input$S
    p_D = input$p_D
    w = as.numeric(strsplit(input$text_w,";")[[1]])
    beta = input$beta
    delta = input$delta
    sigma = input$sigma  

    m = as.numeric(strsplit(input$text_m,";")[[1]])
    
    Z = rmultinom(1,nreads,w)
    reads = experiment(Z,S,m,delta,sigma,beta,p_D)
    reads = GRanges(seqnames = "sim",
      ranges = IRanges(start = reads[,(position)],end = reads[,(position)]),
      strand = ifelse(reads[,(strand)]==1,"+","-"))  
  })


  # make rendering parts

  output$printReads <- renderPrint({
    show(resize(reads(),input$fl))
  })

  output$peak <- renderPlot({
    delta = input$delta
    binding = as.numeric(strsplit(input$text_m,";")[[1]])
    rr = resize(reads(),input$fl)
    fwd_reads = coverage(subset(rr,strand == "+"))[[1]]
    bwd_reads = coverage(subset(rr,strand == "-"))[[1]]
    all = coverage(rr)[[1]]
    fwd = step_fun(fwd_reads,input$S)
    bwd = step_fun(bwd_reads,input$S)
    all = step_fun(all,input$S)
    dt1 = data.table(position = 1:input$S,counts =fwd,strand = "+",panel="strand")
    dt2 = data.table(position = 1:input$S,counts =bwd,strand = "-",panel="strand")
    dt3 = data.table(position = 1:input$S,counts =all,strand = "both",panel="join")
    dt = rbind(dt1,dt2,dt3) 
    dt[,counts:= 1e3*counts / length(rr)]
    p1 = ggplot(dt,aes(position,counts,colour = strand))+
      geom_line(size=1)+
      scale_color_manual(values = c("blue","red","black"))+
      theme(legend.position = "none",axis.text = element_text(size = 12),
            axis.title = element_text(size = 18),
            legend.text = element_text(size =12),
            legend.title = element_text(size= 18),
            strip.text = element_text(size = 18))+
      geom_abline(slope=0,intercept =0,linetype =2)+
      ylab("normalized counts")+    
      geom_vline(xintercept = binding,linetype = 2)+
      geom_vline(xintercept = c(binding - delta,binding + delta),
        linetype = 2,colour = "magenta")+xlab("Genomic coordinates")+
      theme(legend.position = "top")+facet_grid(panel ~. )
    print(p1)
  },height = 700,width = 800)

  output$peakMat <- renderPlot({

    input$simulate

    isolate({
      fr = 25:300
      mat = do.call(rbind,mclapply(fr,fragment_length_cover,reads(),input$S,mc.cores = mc))

      m = as.numeric(strsplit(input$text_m,";")[[1]])     
      x = input$fl
      
      p2 = ggplot(mat,aes(position,fragLen,fill=counts))+geom_tile()+
        scale_fill_gradientn(colours = r)+
        theme(legend.position = "bottom",axis.text = element_text(size = 12),
              axis.title = element_text(size = 18),
              legend.text = element_text(size =12),
              legend.title = element_text(size= 18)
              )+        
        geom_vline(xintercept=m,linetype=2)+ylab("Fragment length")+
        geom_abline(slope=0,intercept = x,colour = "grey",linetype = 3)+
        guides(fill = guide_colorbar(barwidth = unit(400,'points'),
                 barheight = unit(35,'points')))
        print(p2)
    })
          
  },width = 800,height = 680)
  
  ## output$strandRatio <- renderPlot({

  ##   input$calc_fwd

  ##   isolate({
  ##     fr = 25:300
  ##     mat = do.call(rbind,mclapply(fr,forward_strand_ratio,reads(),input$S,mc.cores = mc))

  ##     m = as.numeric(strsplit(input$text_m,";")[[1]])
  ##     x = input$fl

  ##     p2 = ggplot(mat,aes(position,fragLen,fill=fsr))+geom_tile()+
  ##       scale_fill_gradient2(name="Fwd. strand ratio",midpoint=.5)+
  ##       theme(legend.position = "bottom",axis.text = element_text(size = 12),
  ##             axis.title = element_text(size = 18),
  ##             legend.text = element_text(size =12),
  ##             legend.title = element_text(size= 18)
  ##             )+        
  ##       geom_vline(xintercept=m,linetype=2)+ylab("Fragment length")+
  ##       geom_abline(slope=0,intercept = x,colour = "grey",linetype = 3)+
  ##       guides(fill = guide_colorbar(barwidth = unit(400,'points'),
  ##                barheight = unit(35,'points')))
  ##       print(p2)
  ##   })

  ## },width = 800,height = 680)

  
})
