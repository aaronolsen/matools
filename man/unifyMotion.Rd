\name{unifyMotion}
\alias{unifyMotion}
\alias{print.unify_errors}
\title{Unifies marker landmark and motion datasets}
\description{Solves for rigid body transformations for markers fixed in one coordinate system given the motion of those markers in another coordinate system.}
\usage{
unifyMotion(motion, xyz.mat, unify.spec, regexp = FALSE, 
			print.progress = TRUE, print.progress.iter = c(1), verbose = FALSE,
			replace.xyz = TRUE, plot.diag = NULL, cp.use = TRUE, cp.axis.with.vp = FALSE)

\method{print}{unify_errors}(x, ...)
}
\arguments{
  \item{motion}{}
  \item{xyz.mat}{}
  \item{unify.spec}{}
  \item{regexp}{}
  \item{print.progress}{}
  \item{print.progress.iter}{}
  \item{verbose}{}
  \item{replace.xyz}{}
  \item{plot.diag}{}
  \item{cp.use}{}
  \item{cp.axis.with.vp}{}
}
\details{
Details about unifyMotion...
}
\value{
a list containing a motion object  with the following elements:
	\item{motion}{A motion object containing the input marker motion and the calculated rigid body transformations.}
	\item{error}{Unification errors (a list of class \code{"unify_errors"}).}
}
\author{Aaron Olsen}