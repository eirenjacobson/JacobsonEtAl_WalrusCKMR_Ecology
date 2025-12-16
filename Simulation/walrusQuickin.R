#' Look up shared ancestors for `quickin`-style kinship reporting
#'
#' Targets one specific shared ancestor-type, matrix-wise between pairs of
#' animals.
#'
#' @param anc1s a matrix of ancestors for animal 1
#' @param anc2s a matrix of ancestors for animal 2
#' @param same `TRUE`/`FALSE`. `TRUE` if you're looking for shared ancestors at
#' the same position from each member of the pair (e.g., if you're looking for
#' ancestors that are the mother of both pair members.
#' @param weed if non-NULL, will attempt to remove pairs given as the parameter
#' value.

"xpairs" <-
  function( anc1s, anc2s, same, weed=NULL) {
    ## Targeting one specific shared ancestor-type--- eg Mum of 1st == Grandma of 2nd
    ## anc1s & anc2s should be matrices
    # same=TRUE if anc1s == anc2s ie looking for shared mums
    ## weed: any stronger kin to remove from inbred pairs (eg to exclude known POPs from HTPs)--- result of 'both()'
    
    if( FALSE && too_smart) {
      crosses <- c( outer( 1:ncol( anc1s), 1:ncol( anc2s), sprintf( fmt='_%i_%i_')))
      ranc1s <- ...
    }
    
    if( isNamespaceLoaded( 'fastmatch')) {
      match <- fastmatch::fmatch
    }
    
    if( is.null( dim( anc1s))) {
      anc1s <- matrix( anc1s, ncol=1)
    }
    
    if( is.null( dim( anc2s))) {
      anc2s <- matrix( anc2s, ncol=1)
    }
    
    # Find all ancestors that fit the bill
    N <- nrow( anc1s)
    stopifnot( nrow( anc2s)==N)
    
    shareds <- integer()
    for( i1 in 1 %upto% ncol( anc1s)) {
      for( i2 in 1 %upto% ncol( anc2s)) {
        if( same) { # check to see if anc1 matches any anc2 *except* itself
          mm <- match( anc1s[,i1], rev( anc2s[,i2]), 0L)
          mm[ mm != 0] <- N+1-mm[ mm != 0]
          mm[ mm==seq_along( mm)] <- 0
          new_shareds <- anc1s[ mm>0, i1] # should be same as:
          # test <- anc2s[ mm[ mm>0], i1]
        } else {
          combo <- c( anc1s[,i1], anc2s[,i2])
          new_shareds <- combo[ duplicated( combo)]
        }
        
        shareds <- c( shareds, new_shareds)
      }
    }
    # desc1 <- match( anc1s, shareds, 0)
    # desc2 <- match( anc2s, shareds, 0)
    
    # Find all descendent-pairs of each ancestor
    keeps1 <- keeps2 <- integer()
    for( i in shareds) {
      desc1 <- which( rowSums( anc1s==i)>0)
      desc2 <- which( rowSums( anc2s==i)>0)
      keeps1 <- c( keeps1, rep( desc1, length( desc2)))
      keeps2 <- c( keeps2, rep( desc2, each=length( desc1)))
    }
    
    # Zap dups (keep A/B but not B/A or A/A)
    swap <- keeps1 > keeps2
    temp <- keeps1[ swap]
    keeps1[ swap] <- keeps2[ swap]
    keeps2[ swap] <- temp
    
    # now they're in row-wise order, so just check uniqueness of "joint row"--- RHS expression is unique per integer-pair
    both_keep <- both( keeps1, keeps2)
    keep <- (keeps1 != keeps2) & !duplicated( both_keep)
    keeps1 <- keeps1[ keep]
    keeps2 <- keeps2[ keep]
    
    keepair <- cbind( keeps1, keeps2)
    if( !is.null( weed)) {
      keepair <- keepair[ both_keep[ keep] %not.in% weed,,drop=FALSE]
    }
    return( keepair)
  }

#' Unique reals from ordered pairs of integers
#'
#' Return a vector of unique real-valued numbers from a pair of integers
#' @param m1 a vector of integers
#' @param m2 another vector of integers

"both" <-
  function( m1, m2) {
    ## Unique reals from ordered pairs of ints
    if( missing( m2)) {
      m1[,1] + 0.5/m1[,2]
    } else
      m1 + 0.5/m2
  }


#' Quick lookup of CKMR-relevant relationships
#'
#' `quickin` performs quick lookup of the kinships directly relevant to
#' close-kin mark-recapture. It returns a list of eight character
#' arrays, with each array holding one kinship in one pair of animals
#' per row.
#'
#' The named relationship classes (in list order) are:
#' - POPs: Parent-offspring pairs. One is the other's parent.
#' - HSPs: Half-sibling pairs. The individuals share one parent.
#' - FSPs: Full-sibling pairs. The individuals share two parents.
#' - GGPs: Grandparent-grandoffspring pairs. One is the other's grandparent.
#' - HTPs: Half-thiatic pairs. One individual's grandparent is the other
#' individual's parent.
#' - FTPs: Full-thiatic pairs. Two of one individual's grandparents are the
#' other individual's parents.
#' - HCPs: Half-cousin pairs. The individuals share one grandparent.
#' - FCPs: Full-cousin pairs. The individuals share two grandparents.
#'
#' @param inds an `indiv` matrix, as from [mort()], with some individuals
#' marked as 'captured'
#' @param max_gen the maximum depth to look up relatives, in generations.
#' `max_gen = 2` is sufficient for relatives used in CKMR
#' @seealso [fishSim::findRelatives()]
#' @seealso [fishSim::capture()]
#' @export

walrusQuickin <- function (inds, max_gen = 2, cousins. = TRUE, thiatics. = TRUE) 
{
  if (inds %is.not.a% "data.frame") {
    inds <- dfify(inds)
  }
  gc()
  if (isNamespaceLoaded("fastmatch")) {
    match <- fastmatch::fmatch
  }
  keep_me <- with(inds, {
    keep_me <- which(!is.na(SampY))
    prev_keep_me <- keep_me
    for (gen in 1 %upto% max_gen) {
      keep_me_too1 <- match(Mum[prev_keep_me], Me, 0L) %except% 
        keep_me
      keep_me_too2 <- match(Dad[prev_keep_me], Me, 0L) %except% 
        c(keep_me, keep_me_too1)
      prev_keep_me <- c(keep_me_too1, keep_me_too2)
      keep_me <- c(keep_me, prev_keep_me %except% 0L)
    }
    return(unique(keep_me))
  })
  HFsep <- function(H1, H2, keep_sex = FALSE) {
    iFull <- match(both(H1), both(H2), 0)
    Full <- H1[iFull > 0, , drop = FALSE]
    if (nrow(Full)) {
      H1 <- H1[iFull == 0, , drop = FALSE]
      H2 <- H2[-(iFull %except% 0), , drop = FALSE]
    }
    Half <- rbind(H1, H2)
    if (!keep_sex) {
      returnList(Half, Full)
    }
    else {
      returnList(Half, Full, H1, H2)
    }
  }
  extract.named(inds[keep_me, ])
  Sample <- which(!is.na(SampY))
  MOP <- xpairs(Me[Sample], Mum[Sample], same = FALSE)
  FOP <- xpairs(Me[Sample], Dad[Sample], same = FALSE)
  POP <- rbind(MOP, FOP)
  both_so_far <- both(POP)
  MatSP <- xpairs(Mum[Sample], Mum[Sample], same = TRUE, weed = both_so_far)
  PatSP <- xpairs(Dad[Sample], Dad[Sample], same = TRUE, weed = both_so_far)
  {
    HSP
    FSP
    MatHSP
    PatHSP
  } %<-% HFsep(MatSP, PatSP, TRUE)
  isampMum <- match(Mum[Sample], Me)
  isampDad <- match(Dad[Sample], Me)
  GmM <- Mum[isampMum]
  GfM <- Mum[isampDad]
  GM_Sample <- cbind(GmM, GfM)
  GmF <- Dad[isampMum]
  GfF <- Dad[isampDad]
  GF_Sample <- cbind(GmF, GfF)
  both_so_far <- c(both(POP), both(HSP), both(FSP))
  MatGGP <- xpairs(Me[Sample], GM_Sample, same = FALSE, weed = both_so_far)
  PatGGP <- xpairs(Me[Sample], GF_Sample, same = FALSE, weed = both_so_far)
  GGP <- rbind(MatGGP, PatGGP)
  onlyMatGGP <- xpairs(Me[Sample], GM_Sample[, 1], FALSE)
  if (thiatics.) {
    both_so_far <- c(both_so_far, both(GGP))
    MatHTP <- xpairs(Mum[Sample], GM_Sample, FALSE, weed = both_so_far)
    PatHTP <- xpairs(Dad[Sample], GF_Sample, FALSE, weed = both_so_far)
    {
      HTP
      FTP
    } %<-% HFsep(MatHTP, PatHTP, FALSE)
  }
  if (cousins.) {
    both_so_far <- c(both_so_far, both(HTP), both(FTP))
    MatCP <- xpairs(GM_Sample, GM_Sample, same = TRUE, weed = both_so_far)
    PatCP <- xpairs(GF_Sample, GF_Sample, same = TRUE, weed = both_so_far)
    {
      HCP
      FCP
    } %<-% HFsep(MatCP, PatCP, FALSE)
    onlyMatCP <- xpairs(GM_Sample[, 1], GM_Sample[, 1], same = TRUE)
    onlyMatHCP <- HFsep(onlyMatCP, PatCP, TRUE)[[3]]
  }
  ID <- function(m) {
    m[] <- Me[Sample[c(m)]]
    m
  }
  returnees <- cq(MOP, MatHSP)
  returnees <- returnees %that.are.in% lsall(environment())
  retlist <- list()
  for (returnee in returnees) {
    retlist[[returnee]] <- ID(get(returnee))
  }
  return(retlist)
}

