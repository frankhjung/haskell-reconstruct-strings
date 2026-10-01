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

import           Assembler.Internal (assemble)
