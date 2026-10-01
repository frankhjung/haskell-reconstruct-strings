-- |
-- Module      : Assembler.Internal
-- Description : Functional greedy overlap sequence assembler internal implementation
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Pure functional implementation of greedy overlap sequence assembly for
-- reconstructing DNA strands and text from source fragments.
-- This module exposes internal functions for exhaustive unit testing.
module Assembler.Internal
  ( assemble
  , calculateOverlap
  , mergePair
  , findBestOverlap
  , selectBetter
  , compareCandidates
  , filterContainedFragments
  , sortCanonical
  , isProperSubstringOf
  ) where

import           Assembler.Types           (AssemblyError (..), Contig (..),
                                            Fragment (..),
                                            OverlapCandidate (..),
                                            OverlapLength)
import           Data.Containers.ListUtils (nubOrd)
import           Data.List                 (foldl', sortBy)
import qualified Data.Text                 as T

-- | Reconstruct contigs from a collection of fragments given a minimum overlap.
--
-- Returns 'Left' if preconditions are violated ('InvalidMinOverlap' or
-- 'EmptyFragmentEncountered'), otherwise returns 'Right' with contigs sorted
-- canonically.
--
-- >>> assemble [Fragment "ABC", Fragment "BCD", Fragment "CDE"] 2
-- Right [Contig {unContig = "ABCDE"}]
assemble :: [Fragment] -> Int -> Either AssemblyError [Contig]
assemble inputFragments minOverlap
  | minOverlap < 1 = Left (InvalidMinOverlap minOverlap)
  | any (T.null . unFragment) inputFragments = Left EmptyFragmentEncountered
  | otherwise =
      let initialPool   = filterContainedFragments inputFragments
          finalPool     = reducePool initialPool minOverlap
          sortedContigs = sortCanonical (map toContig finalPool)
      in Right sortedContigs

-- | Recursively reduce the candidate pool by greedily merging the best
-- overlap pair and eliminating dynamically contained fragments.
reducePool :: [Fragment] -> Int -> [Fragment]
reducePool [] _ = []
reducePool [sole] _ = [sole]
reducePool pool minOverlap =
      case findBestOverlap pool minOverlap of
        Nothing -> pool
        Just candidate ->
          let merged      = mergePair (prefixFragment candidate)
                                      (suffixFragment candidate)
                                      (matchLength candidate)
              remaining   = [ f | f <- pool
                                , f /= prefixFragment candidate
                                , f /= suffixFragment candidate ]
              updatedPool = filterContainedFragments (merged : remaining)
          in reducePool updatedPool minOverlap

-- | Calculate the longest suffix-prefix overlap between two distinct fragments.
--
-- Returns 0 if no match meeting or exceeding @minOverlap@ is found, or if
-- the match length equals or exceeds the length of the longer fragment.
--
-- >>> calculateOverlap (Fragment "ATGGC") (Fragment "GGCGT") 2
-- 3
calculateOverlap :: Fragment -> Fragment -> Int -> OverlapLength
calculateOverlap (Fragment prefix) (Fragment suffix) minOverlap =
  let maxPossible = min (T.length prefix) (T.length suffix)
      candidates  =
        [ len
        | len <- [maxPossible, maxPossible - 1 .. minOverlap] -- lengths descending
        , T.takeEnd len prefix == T.take len suffix -- test end of prefix with start of suffix
        , len < max (T.length prefix) (T.length suffix) -- proper overlap; not identical
        ]
  in case candidates of
       (best : _) -> best
       []         -> 0

-- | Merge two fragments along an overlapping boundary.
--
-- >>> mergePair (Fragment "ATGG") (Fragment "GGCC") 2
-- Fragment {unFragment = "ATGGCC"}
mergePair :: Fragment -> Fragment -> OverlapLength -> Fragment
mergePair (Fragment prefix) (Fragment suffix) overlapLen =
  Fragment (prefix <> T.drop overlapLen suffix)

-- | Select the preferred overlap candidate according to the assembly ordering.
--
-- The current candidate wins unless the next candidate is strictly better under
-- the deterministic ordering defined by 'compareCandidates'.
--
-- >>> let left = OverlapCandidate (Fragment "AAA") (Fragment "BBB") 2
-- >>> let right = OverlapCandidate (Fragment "CCC") (Fragment "DDD") 3
-- >>> selectBetter left right
-- OverlapCandidate {prefixFragment = Fragment {unFragment = "CCC"}, suffixFragment = Fragment {unFragment = "DDD"}, matchLength = 3}
selectBetter :: OverlapCandidate -> OverlapCandidate -> OverlapCandidate
selectBetter curr nextCandidate =
  case compareCandidates curr nextCandidate of
    LT -> nextCandidate
    _  -> curr

-- | Compare two overlap candidates using the assembly selection rules.
--
-- Ordering is determined by:
-- 1. longer overlap match length first
-- 2. lexicographically smaller prefix fragment first
-- 3. lexicographically smaller suffix fragment first
--
-- >>> let left = OverlapCandidate (Fragment "ABC") (Fragment "XYZ") 3
-- >>> let right = OverlapCandidate (Fragment "ABD") (Fragment "UVW") 3
-- >>> compareCandidates left right
-- GT
compareCandidates :: OverlapCandidate -> OverlapCandidate -> Ordering
compareCandidates a b =
  compare (matchLength a) (matchLength b)
    <> compare (prefixFragment b) (prefixFragment a)
    <> compare (suffixFragment b) (suffixFragment a)

-- | Find the single best overlap candidate across all ordered fragment pairs.
--
-- Implements strict three-tier deterministic tie-breaking:
-- 1. Longest overlap match length (descending)
-- 2. Lexicographically smaller prefix fragment (ascending)
-- 3. Lexicographically smaller suffix fragment (ascending)
--
-- >>> findBestOverlap [Fragment "ABC", Fragment "BCD", Fragment "CDE"] 2
-- Just (OverlapCandidate {prefixFragment = Fragment {unFragment = "ABC"}, suffixFragment = Fragment {unFragment = "BCD"}, matchLength = 2})
findBestOverlap :: [Fragment] -> Int -> Maybe OverlapCandidate
findBestOverlap pool minOverlap =
  let candidates =
        [ OverlapCandidate
            { prefixFragment  = prefix
            , suffixFragment  = suffix
            , matchLength = len
            }
        | prefix <- pool
        , suffix <- pool
        , prefix /= suffix
        , let len = calculateOverlap prefix suffix minOverlap
        , len >= minOverlap
        ]
  in case candidates of
       [] -> Nothing
       (firstCandidate : rest) ->
         Just (foldl' selectBetter firstCandidate rest)

-- | Eliminate exact duplicates and any fragments fully contained as proper
-- substrings inside longer fragments.
--
-- >>> filterContainedFragments [Fragment "ACGT", Fragment "CGT", Fragment "ACGT"]
-- [Fragment {unFragment = "ACGT"}]
filterContainedFragments :: [Fragment] -> [Fragment]
filterContainedFragments fragments =
  filter isNotContained uniqueFragments
  where
    uniqueFragments = nubOrd fragments
    isNotContained f = not (any (isProperSubstringOf f) uniqueFragments)

-- | Check if one fragment is a proper substring of another.
-- Evaluates to True if they are not identical and the first is fully
-- contained in the second.
--
-- >>> isProperSubstringOf (Fragment "CGT") (Fragment "ACGTA")
-- True
isProperSubstringOf :: Fragment -> Fragment -> Bool
isProperSubstringOf s1 s2 =
  s1 /= s2 && unFragment s1 `T.isInfixOf` unFragment s2

-- | Sort contigs into canonical output order:
-- 1. Descending sequence length (longer first)
-- 2. Ascending lexicographical sequence order
--
-- >>> sortCanonical [Contig "DEF", Contig "ABC", Contig "AB"]
-- [Contig {unContig = "DEF"},Contig {unContig = "ABC"},Contig {unContig = "AB"}]
sortCanonical :: [Contig] -> [Contig]
sortCanonical = sortBy compareContigs
  where
    compareContigs (Contig a) (Contig b) =
      compare (T.length b) (T.length a) <> compare a b

-- | Convert an assembled 'Fragment' to a 'Contig'.
toContig :: Fragment -> Contig
toContig = Contig . unFragment
