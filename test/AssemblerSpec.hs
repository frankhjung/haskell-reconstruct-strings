{-# OPTIONS_GHC -Wno-orphans #-}
-- |
-- Module      : AssemblerSpec
-- Description : Unit and property tests for greedy overlap assembler
-- Copyright   : (c) Frank H Jung, 2026
-- License     : BSD-3-Clause
-- Maintainer  : frankhjung@linux.com
-- Stability   : experimental
--
-- Comprehensive test suite validating specification examples, edge cases,
-- tie-breaking rules, containment filtering, and canonical ordering.
module AssemblerSpec (spec) where

import           Assembler       (assemble, calculateOverlap,
                                  filterContainedReads, mergePair)
import           Assembler.Types (AssemblyError (..), Contig (..), Read (..))
import qualified Data.Text       as T
import           Prelude         hiding (Read)
import           Test.Hspec      (Spec, describe, it, shouldBe)
import           Test.QuickCheck (Arbitrary (..), elements, listOf, property)

instance Arbitrary Read where
  arbitrary = Read . T.pack <$> listOf (elements "ACGT")

spec :: Spec
spec = do
  describe "calculateOverlap" $ do
    it "calculates correct overlap length when match exceeds threshold" $
      calculateOverlap (Read "ATGGC") (Read "GGCGT") 2 `shouldBe` 3

    it "returns 0 when overlap is strictly below threshold" $
      calculateOverlap (Read "ATGGC") (Read "CGTGCA") 2 `shouldBe` 0

    it "disallows total containment overlap where match == max length" $
      calculateOverlap (Read "ABC") (Read "ABC") 2 `shouldBe` 0

  describe "mergePair" $
    it "concatenates prefix read with remaining suffix" $
      mergePair (Read "ATGGC") (Read "GGCGT") 3 `shouldBe` Read "ATGGCGT"

  describe "filterContainedReads" $ do
    it "deduplicates identical reads" $
      filterContainedReads [Read "ACGT", Read "ACGT"]
        `shouldBe` [Read "ACGT"]

    it "removes reads that are proper substrings of others" $
      filterContainedReads [Read "ACGT", Read "CGT"]
        `shouldBe` [Read "ACGT"]

  describe "assemble" $ do
    it "assembles single valid overlap (Example A)" $ do
      let inputReads = [Read "ABC", Read "BCD", Read "CDE"]
      assemble inputReads 2 `shouldBe` Right [Contig "ABCDE"]

    it "handles no valid overlap (Example B)" $ do
      let inputReads = [Read "ABC", Read "DEF"]
      assemble inputReads 2
        `shouldBe` Right [Contig "ABC", Contig "DEF"]

    it "removes contained read during pre-processing (Example C)" $ do
      let inputReads = [Read "ACGT", Read "CGT"]
      assemble inputReads 2 `shouldBe` Right [Contig "ACGT"]

    it "eliminates dynamic containment after merge (Example D)" $ do
      let inputReads = [Read "AAATTT", Read "TTTGGG", Read "ATT"]
      assemble inputReads 3 `shouldBe` Right [Contig "AAATTTGGG"]

    it "guarantees canonical permutation invariance (Example E)" $
      assemble [Read "DEF", Read "ABC"] 2
        `shouldBe` assemble [Read "ABC", Read "DEF"] 2

    it "returns error on invalid min overlap" $
      assemble [Read "ACGT"] 0
        `shouldBe` Left (InvalidMinOverlap 0)

    it "returns error on empty read input" $
      assemble [Read "ACGT", Read ""] 2
        `shouldBe` Left EmptyReadEncountered

    it "returns empty contig list on empty input" $
      assemble [] 2 `shouldBe` Right []

    it "returns singleton contig on single read input" $
      assemble [Read "ACGT"] 2 `shouldBe` Right [Contig "ACGT"]

    it "preserves total function safety over arbitrary thresholds" $
      property $ \k ->
        if k < 1
          then assemble [Read "ACGT"] k == Left (InvalidMinOverlap k)
          else assemble [Read "ACGT"] k == Right [Contig "ACGT"]

    it "is permutation invariant (reverse)" $
      property $ \rs k ->
        let k' = max 1 k
         in assemble (reverse rs) k' == assemble rs k'

    it "exhibits idempotence of containment filtering" $
      property $ \rs ->
        filterContainedReads (filterContainedReads rs) == filterContainedReads rs
