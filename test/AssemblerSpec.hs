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

import           Assembler.CLI      (formatAssemblyError, runPipeline,
                                     sanitizeInput)
import           Assembler.Internal (assemble, calculateOverlap,
                                     filterContainedFragments,
                                     isProperSubstringOf, mergePair)
import           Assembler.Types    (AssemblyError (..), Contig (..),
                                     Fragment (..), OverlapCandidate (..))
import qualified Data.Text          as T
import           Test.Hspec         (Spec, describe, it, shouldBe)
import           Test.QuickCheck    (Arbitrary (..), elements, listOf, property)

instance Arbitrary Fragment where
  arbitrary = Fragment . T.pack <$> listOf (elements "ACGT")

instance Arbitrary OverlapCandidate where
  arbitrary = OverlapCandidate <$> arbitrary <*> arbitrary <*> arbitrary

spec :: Spec
spec = do
  describe "calculateOverlap" $ do
    it "calculates correct overlap length when match exceeds threshold" $
      calculateOverlap (Fragment "ATGGC") (Fragment "GGCGT") 2 `shouldBe` 3

    it "returns 0 when overlap is strictly below threshold" $
      calculateOverlap (Fragment "ATGGC") (Fragment "CGTGCA") 2 `shouldBe` 0

    it "disallows total containment overlap where match == max length" $
      calculateOverlap (Fragment "ABC") (Fragment "ABC") 2 `shouldBe` 0

  describe "mergePair" $
    it "concatenates prefix fragment with remaining suffix" $
      mergePair (Fragment "ATGGC") (Fragment "GGCGT") 3 `shouldBe` Fragment "ATGGCGT"

  describe "Ord OverlapCandidate" $ do
    it "prefers the candidate with the longer match length in max" $
      let left  = OverlapCandidate (Fragment "AAA") (Fragment "BBB") 2
          right = OverlapCandidate (Fragment "CCC") (Fragment "DDD") 3
       in max left right `shouldBe` right

    it "breaks ties using the smaller prefix fragment in max" $
      let left  = OverlapCandidate (Fragment "ABC") (Fragment "XYZ") 3
          right = OverlapCandidate (Fragment "ABD") (Fragment "UVW") 3
       in max left right `shouldBe` left

    it "breaks secondary ties using the smaller suffix fragment in max" $
      let left  = OverlapCandidate (Fragment "ABC") (Fragment "UVW") 3
          right = OverlapCandidate (Fragment "ABC") (Fragment "XYZ") 3
       in max left right `shouldBe` left

    it "orders longer overlaps as greater than shorter overlaps" $
      let left  = OverlapCandidate (Fragment "AAA") (Fragment "BBB") 2
          right = OverlapCandidate (Fragment "CCC") (Fragment "DDD") 3
       in compare left right `shouldBe` LT

    it "orders smaller prefix as greater when lengths are equal" $
      let left  = OverlapCandidate (Fragment "ABC") (Fragment "XYZ") 3
          right = OverlapCandidate (Fragment "ABD") (Fragment "UVW") 3
       in compare left right `shouldBe` GT

    it "orders smaller suffix as greater when length and prefix are equal" $
      let left  = OverlapCandidate (Fragment "ABC") (Fragment "UVW") 3
          right = OverlapCandidate (Fragment "ABC") (Fragment "XYZ") 3
       in compare left right `shouldBe` GT

    it "satisfies reflexivity: c <= c" $
      property $ \c -> c <= (c :: OverlapCandidate)

    it "satisfies antisymmetry: a <= b && b <= a ==> a == b" $
      property $ \a b ->
        not (a <= b && b <= (a :: OverlapCandidate)) || a == b

    it "satisfies transitivity: a <= b && b <= c ==> a <= c" $
      property $ \a b c ->
        not (a <= b && b <= (c :: OverlapCandidate)) || a <= c

    it "satisfies totality: a <= b || b <= a" $
      property $ \a b -> a <= b || b <= (a :: OverlapCandidate)

    it "is consistent with Eq: compare a b == EQ <=> a == b" $
      property $ \a b ->
        ((a :: OverlapCandidate) == b) == (a == b)

  describe "filterContainedFragments" $ do
    it "deduplicates identical fragments" $
      filterContainedFragments [Fragment "ACGT", Fragment "ACGT"]
        `shouldBe` [Fragment "ACGT"]

    it "removes fragments that are proper substrings of others" $
      filterContainedFragments [Fragment "ACGT", Fragment "CGT"]
        `shouldBe` [Fragment "ACGT"]

  describe "isProperSubstringOf" $ do
    it "returns True when a fragment is a proper substring of another" $
      isProperSubstringOf (Fragment "CGT") (Fragment "ACGTA") `shouldBe` True

    it "returns False when fragments are identical" $
      isProperSubstringOf (Fragment "ACGT") (Fragment "ACGT") `shouldBe` False

    it "returns False when a fragment is not a substring" $
      isProperSubstringOf (Fragment "CGC") (Fragment "ACGTA") `shouldBe` False

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

    it "exhibits idempotence of containment filtering" $
      property $ \rs ->
        filterContainedFragments (filterContainedFragments rs) == filterContainedFragments rs

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
