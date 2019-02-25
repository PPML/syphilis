\name{initial.files}
\alias{initial.files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
initial.files(state)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{state}{
%%     ~~Describe \code{state} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (state) 
{
    age.dist.dat <<- read.delim("R_inputs/age_nsfg.txt")
    diag.age.sex.dat <<- as.matrix(read.delim(paste("R_inputs/diag_age_sex_rate_", 
        state, ".txt", sep = ""), row.names = 1))
    diag.rr <<- read.delim(paste("R_inputs/case_subpop_rr_", 
        state, ".txt", sep = ""))
    msm.dat <<- as.matrix(read.delim(paste("R_inputs/msm_hiv_proportions_", 
        state, ".txt", sep = ""), row.names = 1))
    stage.dat <<- as.matrix(read.delim(paste("R_inputs/stage_diag_", 
        state, ".txt", sep = ""), row.names = 1))
    p.early.dat <<- read.delim(paste("R_inputs/p_early_", state, 
        ".txt", sep = ""))
    diag.subpop <<- tibble::rownames_to_column(as.data.frame(read.delim(paste("R_inputs/subpop_diag_", 
        state, ".txt", sep = ""), row.names = 1)), "year")
    load(paste("R_inputs/theta_", state, ".rda", sep = ""))
    theta <<- theta
    diag.age.sex.rate.m <<- as.vector(diag.age.sex.dat[, c("diag_all_y_m", 
        "diag_all_o_m")])
    diag.age.sex.rate.m.sd <<- as.vector(diag.age.sex.dat[, c("diag_all_y_m_sd", 
        "diag_all_o_m_sd")])
    diag.age.sex.rate.m.sd <<- diag.age.sex.rate.m.sd * 0.25
    diag.age.sex.rate.f.y <<- as.vector(diag.age.sex.dat[, "diag_all_y_f"])
    diag.age.sex.rate.f.y.sd <<- as.vector(diag.age.sex.dat[, 
        "diag_all_y_f_sd"])
    diag.age.sex.rate.f.o <<- as.vector(diag.age.sex.dat[, "diag_all_o_f"])
    diag.age.sex.rate.f.o.sd <<- as.vector(diag.age.sex.dat[, 
        "diag_all_o_f_sd"])
    rr.diag.subpop <<- as.numeric(diag.rr[, "mean"])
    rr.diag.subpop.sd <<- as.numeric(diag.rr[, "sd"])
    p.msm.cases <<- as.vector(msm.dat[, c("pMSM_y", "pMSM_o")])
    p.msm.cases.var <<- as.vector(msm.dat[, c("pMSM_y_var", "pMSM_o_var")])
    p.hiv.cases <<- as.vector(subset(msm.dat[, c("pHIV_y", "pHIV_o")], 
        (!is.na(msm.dat[, "pHIV_y"])) & (!is.na(msm.dat[, "pHIV_o"]))))
    p.hiv.cases.var <<- as.vector(subset(msm.dat[, c("pHIV_y_var", 
        "pHIV_o_var")], (!is.na(msm.dat[, "pHIV_y_var"])) & (!is.na(msm.dat[, 
        "pHIV_o_var"]))))
    p.sec <<- as.vector(stage.dat[, c("pSec_y_m", "pSec_o_m", 
        "pSec_y_f")])
    p.sec.var <<- as.vector(stage.dat[, c("pSec_var_y_m", "pSec_var_o_m", 
        "pSec_var_y_f")])
    p.el <<- as.vector(stage.dat[, c("pEL_y_m", "pEL_o_m", "pEL_y_f")])
    p.el.var <<- as.vector(stage.dat[, c("pEL_var_y_m", "pEL_var_o_m", 
        "pEL_var_y_f")])
    priors <<- read.delim("R_inputs/priors.txt")
    sd.theta <<- priors$sd.transf.1
    prior.param1 <<- priors$param1
    names(prior.param1) <<- priors$parameter
    prior.param2 <<- priors$param2
    names(prior.param2) <<- priors$parameter
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
