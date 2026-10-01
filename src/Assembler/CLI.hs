-- |
-- Module      : Assembler.CLI
-- Description : Pure transformation pipeline and error formatting for the CLI
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Pure functional pipeline for command-line string reconstruction. Separates
-- text sanitisation, domain error rendering, and sequence assembly from
-- effectful stream I/O.
module Assembler.CLI
  ( sanitizeInput
  , formatAssemblyError
  , runPipeline
  ) where

import           Assembler       (assemble)
import           Assembler.Types (AssemblyError (..), Contig (..),
                                  Fragment (..))
import qualified Data.Text       as T

-- | Sanitise raw input lines by trimming surrounding whitespace and discarding
-- blank or empty lines uniformly across all input sources.
--
-- >>> sanitizeInput ["  ATGGC  ", "", " \t ", "GGCGT"]
-- ["ATGGC","GGCGT"]
sanitizeInput :: [T.Text] -> [T.Text]
sanitizeInput = filter (not . T.null) . map T.strip

-- | Render domain assembly errors into user-facing error text.
--
-- >>> formatAssemblyError (InvalidMinOverlap 0)
-- "Error: minimum overlap must be >= 1 (got 0)"
-- >>> formatAssemblyError EmptyFragmentEncountered
-- "Error: empty fragment encountered in input"
formatAssemblyError :: AssemblyError -> T.Text
formatAssemblyError (InvalidMinOverlap n) =
  "Error: minimum overlap must be >= 1 (got " <> T.pack (show n) <> ")"
formatAssemblyError EmptyFragmentEncountered =
  "Error: empty fragment encountered in input"

-- | Pure end-to-end transformation pipeline from raw text lines and a minimum
-- overlap threshold to either a formatted error message or assembled contig
-- lines.
--
-- >>> runPipeline 2 ["ATGGC", "GGCGT", "CGTGCA"]
-- Right ["ATGGCGTGCA"]
-- >>> runPipeline 0 ["ATGGC"]
-- Left "Error: minimum overlap must be >= 1 (got 0)"
runPipeline :: Int -> [T.Text] -> Either T.Text [T.Text]
runPipeline minOverlap rawLines =
  let cleaned   = sanitizeInput rawLines
      fragments = map Fragment cleaned
  in case assemble fragments minOverlap of
       Left err      -> Left (formatAssemblyError err)
       Right contigs -> Right (map unContig contigs)
