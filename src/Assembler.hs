-- |
-- Module      : Assembler
-- Description : Functional greedy overlap sequence assembler
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Pure functional implementation of greedy overlap sequence assembly for
-- reconstructing DNA strands and text fragments from short reads.
module Assembler
  ( assemble
  , calculateOverlap
  , mergePair
  , findBestOverlap
  , filterContainedReads
  , sortCanonical
  ) where

import           Assembler.Types           (AssemblyError (..), Contig (..),
                                            OverlapCandidate (..),
                                            OverlapLength, Read (..))
import           Data.Containers.ListUtils (nubOrd)
import           Data.List                 (foldl', sortBy)
import qualified Data.Text                 as T
import           Prelude                   hiding (Read, reads)

-- | Reconstruct contigs from a collection of reads given a minimum overlap.
--
-- Returns 'Left' if preconditions are violated ('InvalidMinOverlap' or
-- 'EmptyReadEncountered'), otherwise returns 'Right' with contigs sorted
-- canonically.
--
-- >>> assemble [Read "ABC", Read "BCD", Read "CDE"] 2
-- Right [Contig {unContig = "ABCDE"}]
assemble :: [Read] -> Int -> Either AssemblyError [Contig]
assemble inputReads minOverlap
  | minOverlap < 1 = Left (InvalidMinOverlap minOverlap)
  | any (T.null . unRead) inputReads = Left EmptyReadEncountered
  | otherwise =
      let initialPool   = filterContainedReads inputReads
          finalPool     = reducePool initialPool minOverlap
          sortedContigs = sortCanonical (map toContig finalPool)
      in Right sortedContigs

-- | Recursively reduce the candidate pool by greedily merging the best
-- overlap pair and eliminating dynamically contained reads.
reducePool :: [Read] -> Int -> [Read]
reducePool [] _ = []
reducePool [sole] _ = [sole]
reducePool pool minOverlap =
      case findBestOverlap pool minOverlap of
        Nothing -> pool
        Just candidate ->
          let merged      = mergePair (prefixRead candidate)
                                      (suffixRead candidate)
                                      (matchLength candidate)
              remaining   = [ r | r <- pool
                                , r /= prefixRead candidate
                                , r /= suffixRead candidate ]
              updatedPool = filterContainedReads (merged : remaining)
          in reducePool updatedPool minOverlap

-- | Calculate the longest suffix-prefix overlap between two distinct reads.
--
-- Returns 0 if no match meeting or exceeding @minOverlap@ is found, or if
-- the match length equals or exceeds the length of the longer read.
--
-- >>> calculateOverlap (Read "ATGGC") (Read "GGCGT") 2
-- 3
calculateOverlap :: Read -> Read -> Int -> OverlapLength
calculateOverlap (Read r1) (Read r2) minOverlap =
  let maxPossible = min (T.length r1) (T.length r2)
      candidates  =
        [ len
        | len <- [maxPossible, maxPossible - 1 .. minOverlap] -- desc lengths
        , T.takeEnd len r1 == T.take len r2 -- test end of r1 with start of r2
        , len < max (T.length r1) (T.length r2) -- proper overlap; not identical
        ]
  in case candidates of
       (best : _) -> best
       []         -> 0

-- | Merge two reads along an overlapping boundary.
--
-- >>> mergePair (Read "ATGG") (Read "GGCC") 2
-- Read {unRead = "ATGGCC"}
mergePair :: Read -> Read -> OverlapLength -> Read
mergePair (Read p) (Read s) overlapLen =
  Read (p <> T.drop overlapLen s)

-- | Find the single best overlap candidate across all ordered pairs of reads.
--
-- Implements strict three-tier deterministic tie-breaking:
-- 1. Longest overlap match length (descending)
-- 2. Lexicographically smaller prefix read (ascending)
-- 3. Lexicographically smaller suffix read (ascending)
findBestOverlap :: [Read] -> Int -> Maybe OverlapCandidate
findBestOverlap pool minOverlap =
  let candidates =
        [ OverlapCandidate
            { prefixRead  = r1
            , suffixRead  = r2
            , matchLength = len
            }
        | r1 <- pool
        , r2 <- pool
        , r1 /= r2
        , let len = calculateOverlap r1 r2 minOverlap
        , len >= minOverlap
        ]
  in case candidates of
       [] -> Nothing
       (firstCandidate : rest) ->
         Just (foldl' selectBetter firstCandidate rest)
  where
    selectBetter :: OverlapCandidate -> OverlapCandidate -> OverlapCandidate
    selectBetter curr nextCandidate =
      case compareCandidates curr nextCandidate of
        LT -> nextCandidate
        _  -> curr

    compareCandidates :: OverlapCandidate -> OverlapCandidate -> Ordering
    compareCandidates a b =
      compare (matchLength a) (matchLength b)
        <> compare (prefixRead b) (prefixRead a)
        <> compare (suffixRead b) (suffixRead a)

-- | Eliminate exact duplicates and any reads fully contained as proper
-- substrings inside longer reads.
filterContainedReads :: [Read] -> [Read]
filterContainedReads rawReads =
  let uniqueReads = nubOrd rawReads
  in [ r
     | r <- uniqueReads
     , not (any (\other -> r /= other && unRead r `T.isInfixOf` unRead other)
                uniqueReads)
     ]

-- | Sort contigs into canonical output order:
-- 1. Descending sequence length (longer first)
-- 2. Ascending lexicographical sequence order
sortCanonical :: [Contig] -> [Contig]
sortCanonical = sortBy compareContigs
  where
    compareContigs (Contig a) (Contig b) =
      compare (T.length b) (T.length a) <> compare a b

-- | Convert an assembled 'Read' to a 'Contig'.
toContig :: Read -> Contig
toContig = Contig . unRead
