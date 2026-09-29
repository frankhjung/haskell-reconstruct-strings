-- |
-- Module      : Assembler.Types
-- Description : Domain types for greedy overlap sequence assembly
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Domain models and error types representing genetic fragments, assembled
-- contigs, overlap candidates, and assembly errors.
module Assembler.Types
  ( Fragment (..)
  , Contig (..)
  , OverlapLength
  , OverlapCandidate (..)
  , AssemblyError (..)
  ) where

import           Data.Text (Text)

-- | A single fragment from sequencing data.
newtype Fragment = Fragment
  { unFragment :: Text
  }
  deriving stock (Eq, Ord, Show)

-- | An assembled contiguous sequence resulting from merging fragments.
newtype Contig = Contig
  { unContig :: Text
  }
  deriving stock (Eq, Ord, Show)

-- | Number of characters sharing an exact suffix-prefix match.
type OverlapLength = Int

-- | An ordered candidate pair of fragments with a verified overlap match.
data OverlapCandidate = OverlapCandidate
  { prefixFragment :: !Fragment
    -- ^ The fragment contributing the matching suffix.
  , suffixFragment :: !Fragment
    -- ^ The fragment contributing the matching prefix.
  , matchLength    :: !OverlapLength
    -- ^ Length of the matching suffix-prefix overlap.
  }
  deriving stock (Eq, Ord, Show)

-- | Domain errors returned when precondition validation fails.
data AssemblyError
  = InvalidMinOverlap !Int
    -- ^ The specified minimum overlap threshold is strictly less than 1.
  | EmptyFragmentEncountered
    -- ^ At least one input fragment contains no characters.
  deriving stock (Eq, Show)
