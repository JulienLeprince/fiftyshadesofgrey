library(cstmr)
library(stringr)
library(lubridate)
library(scales)

path_in = "fiftyshadesofgrey/data/out/"
path_out = "fiftyshadesofgrey/fig/"

# Set the working directory. Change this to the location of the example on the computer. Note that "/" is always used in R, also in Windows
setwd("fiftyshadesofgrey/src/")

# Load processed results
df_res <- read.csv( paste(path_in,'all_greybox_fits.csv', sep=""))
df_paths <- read.csv( paste(path_in,'all_greybox_paths.csv', sep=""))


# -------------------------------------------- Calculating CP, CPBE and nCPBES for model evaluation --------------------------------------------
# Get all input files names
list_file_names <- list.files(pattern = "*_finalfit.rda")

# Initialize results lists
x_list <- list()
y_list <- list()
CP_list<- list()
BE_list <- list()
nCPBES_list <- list()

for (file in list_file_names) {
  # Get blg uuid from file name
  split_file_name <- str_split_fixed(file, "_", 3)
  uuid_filtering <- sapply(split_file_name, tail, 1)[[1]]
  # Load results model
  load(file)
  
  # Testing if the result evaluated is correct
  if (length(best_fit) == 1) {
    next
  }

  # Loading data associated with results model
  df <- read.csv(paste(path_in,"blg_",uuid_filtering,".csv",sep=""),sep=";",header=TRUE)

  # Convert to datetime
  df$timedate <- ymd_hms(df$t)
  ## df$t is now hours since start of the experiment.
  df$t <- seq(0, length(df$t)-1, by=1) / 4  # dt = 15 minutes

  # Define input X
  X <- df[c("t", "Ps", "Ta", "Thm", "yTi", "timedate")]
  # Removing NA values
  X = X[complete.cases(X), ]
  ## ----oneStepPred---------------------------------------------------------
  ## Calculate the one-step predictions of the state (i.e. the residuals)
  tmp <- predict(best_fit)[[1]]
  ## Calculate the residuals and put them with the data in a data.frame X
  X$residuals <- X$yTi - tmp$output$pred$yTi
  X$yTiHat <- tmp$output$pred$yTi
  ts <- X$residuals
  # Periodogram calculation
  nCPBES, BE, CP, x, xm, crit <- nCPBE_calc(ts, all=TRUE)
  
  # Saving results to lists
  x_list <- append (x_list, list(x))
  CP_list <- append (CP_list, list(CP))
  BE_list <- append (BE_list, list(BE))
  nCPBES_list <- append (nCPBES_list, list(nCPBE))
}




# --------------------- Plotting obtained Cumulated Periodogram, Boundary Excess and nCPBES ---------------------

# Defining path and picture name for saving
png(paste(path_out,'nCPBES_all.png', sep=""), width=1200, height=480, res=200)
par(mfrow=c(1,3))

# Colors initialization
color_good <- alpha(rgb(0.6196078431372549, 0.00392156862745098, 0.25882352941176473), 0.08)
color_ok <- alpha(rgb(0.998077662437524, 0.9992310649750096, 0.7460207612456747), 0.08)
color_poor <- alpha(rgb(0.3686274509803922, 0.30980392156862746, 0.6352941176470588), 0.08)

### Cumulated Periodogram Plot
for (i in 1:length(list_file_names)) {
  # Color parameter adaptation
  nCPBES <- tail(unlist(nCPBES_list[i], use.names=FALSE), n=1)
  if (is.null(nCPBES)) {
    next
  }else if (nCPBES < 0.01) {
    color <- color_good
  } else if (nCPBES < 0.03) {
    color <- color_ok
  } else {
    color <- color_poor
  }
  
  if (i == 1) {
      fig <- plot(unlist(x_list[i], use.names=FALSE), unlist(CP_list[i], use.names=FALSE), type="s", xlim=c(0, xm),
                  ylim=c(0, 1), xaxs="i", yaxs="i", xlab="Frequency (2/h)",
                  ylab="", col=color)
  } else {
     fig + lines(unlist(x_list[i], use.names=FALSE), unlist(CP_list[i], use.names=FALSE), type="s", xlim=c(0, xm),
                 ylim=c(0, 1), xaxs="i", yaxs="i", xlab="Frequency (2/h)",
                 ylab="", col=color) 
  }
}
lines(c(0, xm*(1-crit)), c(crit, 1), col = "blue", lty = 2)
lines(c(xm*crit, xm), c(0, 1-crit), col = "blue", lty = 2)
title(main = "Cumulated Periodogram")
invisible()

### Boundary Excess Plot
for (i in 1:length(list_file_names)) {
  # Color parameter adaptation
  nCPBES <- tail(unlist(nCPBES_list[i], use.names=FALSE), n=1)
  if (is.null(nCPBES)) {
    next
  }else if (nCPBES < 0.01) {
    color <- color_good
  } else if (nCPBES < 0.03) {
    color <- color_ok
  } else {
    color <- color_poor
  }
  if (i == 1) {
      fig <- plot(unlist(x_list[i], use.names=FALSE), unlist(BE_list[i], use.names=FALSE), type="s", xlim=c(0, xm),
                  ylim=c(0, 1), xaxs="i", yaxs="i", xlab="Frequency (2/h)",
                  ylab="", col=color)
  } else {
     fig + lines(unlist(x_list[i], use.names=FALSE), unlist(BE_list[i], use.names=FALSE), type="s", xlim=c(0, xm),
                 ylim=c(0, 1), xaxs="i", yaxs="i", xlab="Frequency (2/h)",
                 ylab="", col=color) 
  }
}
title(main = "Boundary Excess")
invisible()

### nCPBE Plot
for (i in 1:length(list_file_names)) {
  # Color change
  nCPBES <- tail(unlist(nCPBES_list[i], use.names=FALSE), n=1)
  if (is.null(nCPBES)) {
    next
  }else if (nCPBES < 0.01) {
    color <- color_good
  } else if (nCPBES < 0.03) {
    color <- color_ok
  } else {
    color <- color_poor
  }
  if (i == 1) {
      fig <- plot(unlist(x_list[i], use.names=FALSE), unlist(nCPBES_list[i], use.names=FALSE), type="s", xlim=c(0, xm),
                  ylim=c(0, .1), xaxs="i", yaxs="i", xlab="Frequency (2/h)",
                  ylab="", col=color)
  } else {
     fig + lines(unlist(x_list[i], use.names=FALSE), unlist(nCPBES_list[i], use.names=FALSE), type="s", xlim=c(0, xm),
                 ylim=c(0, .1), xaxs="i", yaxs="i", xlab="Frequency (2/h)",
                 ylab="", col=color) 
  }
}
text(x = 0.37, y = .095, "Poor fits")
lines(c(0, 1), c(.03, .03), col = "red", lty = 2)
text(x = 0.4, y = .025, "Close fits")
lines(c(0, 1), c(.01, .01), col = "red", lty = 2)
text(x = 0.4, y = .005, "Good fits")

title(main = "Cumulated nCPBES")
invisible()
dev.off()


