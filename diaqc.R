#R

## Christian Panse <cp@fgcz.ethz.ch>
## 20231024
## diaqc R

stopifnot(require(readr), require(reshape2))

S <- readr::read_delim("~/Downloads/output.txt", 
                     delim = ";", escape_double = FALSE, col_types = cols(Time = col_datetime(format = "%s"), 
                                                                          Size = col_integer(), Precursors.Identified = col_integer(), 
                                                                          Proteins.Identified = col_integer()), 
                     trim_ws = TRUE) 



S$Instrument <- NA
.assignInstrument <- function(x){
  for (i in c('LUMOS_2', 'EXPLORIS_2')){
   idx <- grepl(i, x$File.Name) 
   x$Instrument[which(idx)] <- i
  }
  x
}
S |> .assignInstrument() -> S

S |> reshape2::melt(id.vars = c("Md5", "File.Name", "Time", "Instrument")) -> SS 

# S |> reshape2::melt(id.vars = c("Md5", "File.Name", "Time", "Size", "Instrument")) -> SS 


pdf("~/Downloads/diannqc.pdf", 5,  30)
lattice::xyplot(value ~ Time | variable,
                group = Instrument,
                data = SS,
                scales = 'free',
                type = 'b',
                layout = c(1,18),
                auto.key = list(space = "bottom"))
dev.off()
