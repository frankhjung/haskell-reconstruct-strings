-- |
-- Module      : Assembler.Types
-- Description : Domain types for greedy overlap sequence assembly
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Domain models and error types representing genetic reads, assembled
-- contigs, overlap candidates, and assembly errors.
module Assembler.Types
  ( Read (..)
  , Contig (..)
  , OverlapLength
  , OverlapCandidate (..)
  , AssemblyError (..)
  ) where

import           Data.Text (Text)
import           Prelude   hiding (Read)

-- | A single sequence fragment read from sequencing data.
newtype Read = Read
  { unRead :: Text
  }
  deriving stock (Eq, Ord, Show)

-- | An assembled contiguous sequence resulting from merging reads.
newtype Contig = Contig
  { unContig :: Text
  }
  deriving stock (Eq, Ord, Show)

-- | Number of characters sharing an exact suffix-prefix match.
type OverlapLength = Int

-- | An ordered candidate pair of reads with a verified overlap match.
data OverlapCandidate = OverlapCandidate
  { prefixRead  :: !Read
    -- ^ The read contributing the matching suffix.
  , suffixRead  :: !Read
    -- ^ The read contributing the matching prefix.
  , matchLength :: !OverlapLength
    -- ^ Length of the matching suffix-prefix overlap.
  }
  deriving stock (Eq, Ord, Show)

-- | Domain errors returned when precondition validation fails.
data AssemblyError
  = InvalidMinOverlap !Int
    -- ^ The specified minimum overlap threshold is strictly less than 1.
  | EmptyReadEncountered
    -- ^ At least one input read contains no characters.
  deriving stock (Eq, Show)
