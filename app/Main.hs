-- |
-- Module      : Main
-- Description : Command-line interface for greedy overlap sequence assembler
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Executable entry point accepting sequence fragments and assembling them
-- into canonical contigs using greedy overlap reduction.
module Main (main) where

import           Assembler.CLI       (runPipeline)
import qualified Data.Text           as T
import qualified Data.Text.IO        as TIO
import           Options.Applicative (Parser, auto, execParser, fullDesc,
                                      header, help, helper, info, long, metavar,
                                      option, optional, progDesc, short,
                                      showDefault, strOption, value)
import           System.Exit         (exitFailure)
import           System.IO           (stderr)

-- | Command-line configuration parameters.
data Options = Options
  { optMinOverlap :: !Int
  , optInputFile  :: !(Maybe FilePath)
  }

optionsParser :: Parser Options
optionsParser =
  Options
    <$> option
          auto
          ( long "min-overlap"
              <> short 'm'
              <> metavar "INT"
              <> value 2
              <> showDefault
              <> help "Minimum overlap threshold"
          )
    <*> optional
          ( strOption
              ( long "file"
                  <> short 'f'
                  <> metavar "FILE"
                  <> help "Input file containing one fragment per line"
              )
          )

-- | Run CLI application.
main :: IO ()
main = do
  opts <- execParser optsInfo
  rawLines <- maybe (T.lines <$> TIO.getContents) -- default read from stdin
                    (fmap T.lines . TIO.readFile) -- from a file
                    (optInputFile opts)           -- maybe a file name
  either reportError
         (mapM_ TIO.putStrLn)
         (runPipeline (optMinOverlap opts) rawLines)
  where
    reportError err = TIO.hPutStrLn stderr err >> exitFailure
    optsInfo =
      info
        (helper <*> optionsParser)
        ( fullDesc
            <> progDesc
              "Reconstruct DNA sequences from overlapping fragment sequences"
            <> header "reconstruct-strings - greedy overlap sequence assembler"
        )
