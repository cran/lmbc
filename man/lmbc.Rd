\name{lmbc}
\alias{lmbc}
\title{Linear Model Bias Correction for RNA-Seq Data}
\description{This function implements the method in "Bias Correction in RNA-Seq Short-read Counts using Penalized Regression." It first using a penalized regression to determine an appropriate surrounding sequence, then refits the model using the dinucleotide expansion.}
\usage{lmbc(data, up, down, power)}
\arguments{
	\item{data}{Data in the format accepted by the mseq package, for use with its expansion function.}
	\item{up}{Length of the initial upstream sequence.}
	\item{down}{Length of the initial downstream sequence.}
	\item{power}{Additional parameter for weighted fit penalty. Default and recommended is 3.}
}

\value{
	\item{seqLen}{Length of resulting sequence.}
	\item{lasL}{Length of resulting upstream sequence.}
	\item{lasR}{Length of resulting downstream sequence.}
	\item{r2}{Coefficient of determination for the model.}
	\item{logLik}{Log-likelihood of the model.}
}

\author{David Dalpiaz}


\references{
    Dalpiaz, D., He, X., and Ma, P. (2012) Bias correction in RNA-Seq short-read counts using penalized regression , Statistics in Biosciences , DOI: 10.1007/s12561-012-9057-6.
}

\examples{
    data(g1_part)
    lmbc.ex <- lmbc(g1_part,40,41,3)
    lmbc.ex
}
