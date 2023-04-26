library(Cardinal)
library(tidyr)
library(ggplot2)
library(gridExtra)

select_cells_from_single_mz <- function(f, mz,i) {
  # f is the MSI data read by readMSIData
  # mz is the m/z that can be used to distinguish cells
  # i is the intensity cutoff that can be used for filtering
  df <- as.data.frame(as.matrix(spectra(f)))
  XCoord <- coord(f)$x
  YCoord <- coord(f)$y
  coord <- data.frame(x=XCoord,y=YCoord)
  df <- df[features(f, mz=mz),]
  df <- df %>%
    pivot_longer(
      cols = starts_with("V"),
      values_to = "intensity",
    )
  df <- df[,2]
  df <- cbind(coord,df)
  df2 <- subset(df,intensity>i)
  p1 <- ggplot(df, aes(x, y, fill=intensity)) + geom_tile()
  p2 <- ggplot(df2, aes(x, y, fill=intensity)) + geom_tile()
  grid.arrange(p1, p2, ncol = 2)
}


f <- readMSIData("./data_analysis/xd_try_01.imzML")
select_cells_from_single_mz(f=f,mz=283.2657,i=800)



