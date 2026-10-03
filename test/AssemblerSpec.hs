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

import           Assembler       (assemble)
import           Assembler.CLI   (formatAssemblyError, runPipeline,
                                  sanitizeInput)
import           Assembler.Types (AssemblyError (..), Contig (..),
                                  Fragment (..))
import qualified Data.Text       as T
import           Test.Hspec      (Spec, describe, it, shouldBe)
import           Test.QuickCheck (Arbitrary (..), elements, listOf, property)

instance Arbitrary Fragment where
  arbitrary = Fragment . T.pack <$> listOf (elements "ACGT")

spec :: Spec
spec = do
  describe "assemble" $ do
    it "assembles single valid overlap (Example A)" $ do
      let inputFragments = [Fragment "ABC", Fragment "BCD", Fragment "CDE"]
      assemble inputFragments 2 `shouldBe` Right [Contig "ABCDE"]

    it "handles no valid overlap (Example B)" $ do
      let inputFragments = [Fragment "ABC", Fragment "DEF"]
      assemble inputFragments 2
        `shouldBe` Right [Contig "ABC", Contig "DEF"]

    it "removes contained fragment during pre-processing (Example C)" $ do
      let inputFragments = [Fragment "ACGT", Fragment "CGT"]
      assemble inputFragments 2 `shouldBe` Right [Contig "ACGT"]

    it "eliminates dynamic containment after merge (Example D)" $ do
      let inputFragments = [Fragment "AAATTT", Fragment "TTTGGG", Fragment "ATT"]
      assemble inputFragments 3 `shouldBe` Right [Contig "AAATTTGGG"]

    it "guarantees canonical permutation invariance (Example E)" $
      assemble [Fragment "DEF", Fragment "ABC"] 2
        `shouldBe` assemble [Fragment "ABC", Fragment "DEF"] 2

    it "returns error on invalid min overlap" $
      assemble [Fragment "ACGT"] 0
        `shouldBe` Left (InvalidMinOverlap 0)

    it "returns error on empty fragment input" $
      assemble [Fragment "ACGT", Fragment ""] 2
        `shouldBe` Left EmptyFragmentEncountered

    it "returns empty contig list on empty input" $
      assemble [] 2 `shouldBe` Right []

    it "returns singleton contig on single fragment input" $
      assemble [Fragment "ACGT"] 2 `shouldBe` Right [Contig "ACGT"]

    it "preserves total function safety over arbitrary thresholds" $
      property $ \k ->
        if k < 1
          then assemble [Fragment "ACGT"] k == Left (InvalidMinOverlap k)
          else assemble [Fragment "ACGT"] k == Right [Contig "ACGT"]

    it "is permutation invariant (reverse)" $
      property $ \rs k ->
        let k' = max 1 k
         in assemble (reverse rs) k' == assemble rs k'



  describe "Assembler.CLI" $ do
    describe "sanitizeInput" $ do
      it "strips surrounding whitespace from lines" $
        sanitizeInput ["  ATGGC  ", "GGCGT\t"] `shouldBe` ["ATGGC", "GGCGT"]

      it "discards empty and whitespace-only lines" $
        sanitizeInput ["", "  ", "\t\n", "ATGGC"] `shouldBe` ["ATGGC"]

    describe "formatAssemblyError" $ do
      it "formats InvalidMinOverlap error message" $
        formatAssemblyError (InvalidMinOverlap 0)
          `shouldBe` "Error: minimum overlap must be >= 1 (got 0)"

      it "formats EmptyFragmentEncountered error message" $
        formatAssemblyError EmptyFragmentEncountered
          `shouldBe` "Error: empty fragment encountered in input"

    describe "runPipeline" $ do
      it "assembles raw lines into formatted contigs" $
        runPipeline 2 ["  ATGGC  ", "", "GGCGT", "CGTGCA"]
          `shouldBe` Right ["ATGGCGTGCA"]

      it "returns formatted error when minOverlap is invalid" $
        runPipeline 0 ["ATGGC"]
          `shouldBe` Left "Error: minimum overlap must be >= 1 (got 0)"

      it "returns empty output for empty input" $
        runPipeline 2 [] `shouldBe` Right []
