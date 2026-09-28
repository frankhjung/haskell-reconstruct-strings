-- |
-- Module      : Main
-- Description : Command-line interface for greedy overlap sequence assembler
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Executable entry point accepting short read fragments and assembling them
-- into canonical contigs using greedy overlap reduction.
module Main (main) where

import           Assembler           (assemble)
import           Assembler.Types     (AssemblyError (..), Contig (..),
                                      Read (..))
import           Control.Applicative (many)
import qualified Data.Text           as T
import qualified Data.Text.IO        as TIO
import           Options.Applicative (Parser, auto, execParser, fullDesc,
                                      header, help, helper, info, long, metavar,
                                      option, optional, progDesc, short,
                                      showDefault, strArgument, strOption,
                                      value)
import           Prelude             hiding (Read)
import           System.Exit         (exitFailure)
import           System.IO           (stderr)

-- | Command-line configuration parameters.
data Options = Options
  { optMinOverlap :: !Int
  , optInputFile  :: !(Maybe FilePath)
  , optReads      :: ![String]
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
                  <> help "Input file containing one read per line"
              )
          )
    <*> many
          ( strArgument
              ( metavar "READ..."
                  <> help "Read fragments passed as arguments"
              )
          )

-- | Run CLI application.
main :: IO ()
main = do
  opts <- execParser optsInfo
  inputContent <- readInput opts
  let readsList = map (Read . T.strip) inputContent
  case assemble readsList (optMinOverlap opts) of
    Left (InvalidMinOverlap n) -> do
      TIO.hPutStrLn stderr $ "Error: minimum overlap must be >= 1 (got "
                          <> T.pack (show n)
                          <> ")"
      exitFailure
    Left EmptyReadEncountered -> do
      TIO.hPutStrLn stderr "Error: empty read encountered in input"
      exitFailure
    Right contigs ->
      mapM_ (TIO.putStrLn . unContig) contigs
  where
    optsInfo =
      info
        (helper <*> optionsParser)
        ( fullDesc
            <> progDesc
              "Reconstruct DNA sequences from overlapping read fragments"
            <> header "reconstruct-strings - greedy overlap sequence assembler"
        )

    readInput :: Options -> IO [T.Text]
    readInput opts = case (optInputFile opts, optReads opts) of
      (Just filePath, _) ->
        filter (not . T.null) . T.lines <$> TIO.readFile filePath
      (Nothing, rs)
        | not (null rs) -> pure (map T.pack rs)
        | otherwise     -> filter (not . T.null) . T.lines <$> TIO.getContents
