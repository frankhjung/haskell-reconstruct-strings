-- |
-- Module      : Assembler
-- Description : Functional greedy overlap sequence assembler
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Pure functional implementation of greedy overlap sequence assembly for
-- reconstructing DNA strands and text from source fragments.
-- This public module exports only the main assembly contract.
module Assembler
  ( assemble
  ) where

import           Assembler.Types           (AssemblyError (..), Contig (..),
                                            Fragment (..), OverlapLength)
import           Data.Containers.ListUtils (nubOrd)
import           Data.List                 (sortBy)
import qualified Data.Text                 as T

-- | An ordered candidate pair of fragments with a verified overlap match.
data OverlapCandidate = OverlapCandidate
  { prefixFragment :: !Fragment
    -- ^ The fragment contributing the matching suffix.
  , suffixFragment :: !Fragment
    -- ^ The fragment contributing the matching prefix.
  , matchLength    :: !OverlapLength
    -- ^ Length of the matching suffix-prefix overlap.
  }
  deriving stock (Eq, Show)

-- | Strict three-tier deterministic total ordering:
-- 1. Longest overlap match length (descending)
-- 2. Lexicographically smaller prefix fragment (ascending)
-- 3. Lexicographically smaller suffix fragment (ascending)
instance Ord OverlapCandidate where
  compare a b =
    compare (matchLength a) (matchLength b)
      <> compare (prefixFragment b) (prefixFragment a)
      <> compare (suffixFragment b) (suffixFragment a)

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
mergePair :: Fragment -> Fragment -> OverlapLength -> Fragment
mergePair (Fragment prefix) (Fragment suffix) overlapLen =
  Fragment (prefix <> T.drop overlapLen suffix)

-- | Find the single best overlap candidate across all ordered fragment pairs.
--
-- Implements strict three-tier deterministic tie-breaking via 'Ord':
-- 1. Longest overlap match length (descending)
-- 2. Lexicographically smaller prefix fragment (ascending)
-- 3. Lexicographically smaller suffix fragment (ascending)
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
       cs -> Just (maximum cs)

-- | Eliminate exact duplicates and any fragments fully contained as proper
-- substrings inside longer fragments.
filterContainedFragments :: [Fragment] -> [Fragment]
filterContainedFragments fragments =
  filter isNotContained uniqueFragments
  where
    uniqueFragments = nubOrd fragments
    isNotContained f = not (any (isProperSubstringOf f) uniqueFragments)

-- | Check if one fragment is a proper substring of another.
-- Evaluates to True if they are not identical and the first is fully
-- contained in the second.
isProperSubstringOf :: Fragment -> Fragment -> Bool
isProperSubstringOf s1 s2 =
  s1 /= s2 && unFragment s1 `T.isInfixOf` unFragment s2

-- | Sort contigs into canonical output order:
-- 1. Descending sequence length (longer first)
-- 2. Ascending lexicographical sequence order
sortCanonical :: [Contig] -> [Contig]
sortCanonical = sortBy compareContigs
  where
    compareContigs (Contig a) (Contig b) =
      compare (T.length b) (T.length a) <> compare a b

-- | Convert an assembled 'Fragment' to a 'Contig'.
toContig :: Fragment -> Contig
toContig = Contig . unFragment
