background_reads <- function(strand,beta,S)
{
  out = rep(0,length(strand))
  out[strand == 0] = runif(sum(strand == 0),min=1,max=S+beta-1)
  out[strand == 1] = runif(sum(strand == 1),min = 2 -beta, max = S)
  return(round(out,0))
}

chip_reads <- function(strand, m , delta, sigma)
{
  out = rep(0,length(strand))
  out[strand == 0] = rnorm(sum(strand == 0),mean = m + delta,sd =sigma)
  out[strand == 1] = rnorm(sum(strand == 1),mean = m - delta,sd =sigma)
  return(round(out,0))
}

experiment <- function(Z,S,m,delta,sigma,beta,strand_prob)
{
  stopifnot(length(Z) == length(m) + 1)
  strands = lapply(Z,rbinom,1,strand_prob) 
  background = background_reads(strands[[1]],beta,S)
  chip = mapply(chip_reads,strands[-1],m,MoreArgs = list(delta,sigma),SIMPLIFY=FALSE)
  chip = mapply(function(cc,strand,type){
    return(data.table(position = cc,strand = strand,type = paste0("chip",type)))
  },chip,strands[-1],1:length(m),SIMPLIFY=FALSE)
  dt1 = data.table(position = background,strand = strands[[1]],type = "background")
  dt2 = do.call(rbind,chip)
  return(rbind(dt1,dt2))  
}


step_fun <- function(cover,S)
{
  xi = runLength(cover)
  yi = runValue(cover)
  y = stepfun(cumsum(xi),c(0,yi))(1:S)
  return(y)
}


fix_length <- function(x,M)
{
  values = runValue(x)
  lengths = runLength(x)
  if(sum(lengths) < M)
  {
  values = c(values,0)
  lengths = c(lengths,M - sum(lengths))
  }
  out = Rle(values,lengths)
  return(out)  
}


fragment_length_cover<- function(fl,gr,S)
{
  gr = resize(gr,fl)
  all = coverage(gr)[[1]]
  all = step_fun(all,S)
  dt = data.table(position = 1:S,counts =all/fl)
  dt[,fragLen:=fl]
  return(dt)
}
