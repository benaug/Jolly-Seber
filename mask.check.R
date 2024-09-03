mask.check <- function(dSS=NA,cells=NA,n.cells=NA,n.cells.x=NA,n.cells.y=NA,res=NA,xlim=NA,ylim=NA,
                       x.vals=NA,y.vals=NA){
  
  if(nrow(dSS)!=n.cells)stop("'dSS' must have 'n.cells' rows")
  if(ncol(dSS)!=2)stop("'dSS' must have 2 columns")
  if(dim(cells)[1]!=n.cells.x)stop("'cells' should be of length 'n.cells.x'")
  if(dim(cells)[2]!=n.cells.y)stop("'cells' should be of length 'n.cells.y'")
  if(!all(range(x.vals)+c(-res/2,res/2)==xlim))stop("x.vals doesn't match up with xlim")
  if(!all(range(y.vals)+c(-res/2,res/2)==ylim))stop("y.vals doesn't match up with ylim")
  
  for(i in 1:n.cells){
    s.cell.x <- i%%n.cells.x 
    s.cell.y <- floor(i/n.cells.x)+1
    if(s.cell.x==0){
      s.cell.x=n.cells.x
      s.cell.y=s.cell.y-1
    }
    match=which(cells==i,arr.ind=TRUE)
    if(!all(c(s.cell.x,s.cell.y)==match))stop("error 1")
    
    xlim.cell=c(s.cell.x-1,s.cell.x)*res
    ylim.cell=c(s.cell.y-1,s.cell.y)*res
    
    if(max(x.vals[s.cell.x]+c(-res/2,res/2)-xlim.cell)>1e-6)stop("error 2")
    if(max(y.vals[s.cell.y]+c(-res/2,res/2)-ylim.cell)>1e-6)stop("error 3")
  }
  print("all tests passed")
}
