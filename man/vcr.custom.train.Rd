\name{vcr.custom.train}
\alias{vcr.custom.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of a classification made by a custom algorithm on training data.
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{silplot}, \code{qresplot}, \code{mdscolorscale}. The user first needs to train its classfier of choice on the training set. Then here the predicted posteriors probabilities on the training set are needed along with the true labels. If a distToClass matrix is provided, the function calculates 'farness', needed to enable \code{classmap} visualization.
}

\usage{
vcr.custom.train(y, probs, distToClass=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{y}{factor with the given/true class labels of the objects in the training set considered. Make sure that the levels are in the same order as used in training the classifier for a properly order output of the posterior probabilities used in \code{probs}. It's advisable to factorize the response \code{y} before feeding it in the classifier training process. (if another framework, other than R is used, or the algorithm do not accept response factors, take particular attention to the levels ordering, for example order of levels should be the same of the columns of its binary "once-hot-encoded" \code{y} response vector.)}
\item{probs}{Matrix cointaining the posterior probabilities obtained on the considered training set by the trained classifier. Each row represents the set of posterior probabilities for a certain observation for the different classes, that are represented in the columns. The order of the observations in the rows must match the same order used in \code{y}. The classes in the columns must be in the same order as the levels of \code{y}. The posteriors probabilities should be consistent with the fact that the classifier assign the observation to the class for which the posterior probability is higher. Permits 'PAC' computation. Each row should add up to 1.}
\item{distToclass}{Matrix containing the distances of each observation of the training set considered to each possible class. The user should produce this matrix using a custom-devised metric that is coherent with how the trained classifier works and sees the data. The ordering for rows and columns is the same as in \code{probs}. It permits 'farness' computation for each observation and enables the \code{classmap} visual.}

}


\value{
A list with components: \cr
  \item{yint}{number of the given class of each case. Can contain \code{NA}'s.}
  \item{y}{given class label of each case. Can contain \code{NA}'s.}
  \item{levels}{levels of \code{y}}
  \item{predint}{predicted class number of each case. For each case this is the class with the highest posterior probability. Always exists.}
  \item{pred}{predicted label of each case.}
  \item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{y} is missing.}
  \item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{PAC}{probability of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{PCAfits}{if not \code{NULL}, PCA fits to each class, estimated from the training data but also useful for new data.}
  \item{figparams}{parameters for computing \code{fig}, can be used for new data. Is \code{NULL} if distToClass is not feeded. }
  \item{fig}{distance of each case \eqn{i} from each class \eqn{g}. Is \code{NULL} if distToClass is not feeded.}
  \item{farness}{farness of each case from its given class. Is \code{NULL} if distToClass is not feeded.}
  \item{ofarness}{for each case \code{i}, its lowest  \code{fig[i,g]} to any class \code{g}. Is \code{NULL} if distToClass is not feeded.}
}
\references{--add paper at release}
}

\author{
 Lazzari. L
}

\seealso{
\code{\link{vcr.custom.newdata}}, \code{\link{classmap:silplot}}, \code{\link{classmap:qresplot}}, \code{\link{mdscolorscale}}
}

\examples{
# For examples, we refer to the vignette:
\dontrun{
vignette("MAKE EXAMPLE VIGNETTE")
}
}