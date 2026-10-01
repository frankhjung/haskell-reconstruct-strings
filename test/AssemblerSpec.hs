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

import           Assembler.Internal (assemble, calculateOverlap, compareCandidates,
                                     filterContainedFragments, isProperSubstringOf,
                                     mergePair, selectBetter)
import           Assembler.Types (AssemblyError (..), Contig (..),
                                  Fragment (..), OverlapCandidate (..))
import qualified Data.Text       as T
import           Test.Hspec      (Spec, describe, it, shouldBe)
import           Test.QuickCheck (Arbitrary (..), elements, listOf, property)

instance Arbitrary Fragment where
  arbitrary = Fragment . T.pack <$> listOf (elements "ACGT")

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

  describe "selectBetter" $ do
    it "prefers the candidate with the longer match length" $
      let left  = OverlapCandidate (Fragment "AAA") (Fragment "BBB") 2
          right = OverlapCandidate (Fragment "CCC") (Fragment "DDD") 3
       in selectBetter left right `shouldBe` right

    it "breaks ties using the smallest prefix fragment" $
      let left  = OverlapCandidate (Fragment "ABC") (Fragment "XYZ") 3
          right = OverlapCandidate (Fragment "ABD") (Fragment "UVW") 3
       in selectBetter left right `shouldBe` left

  describe "compareCandidates" $ do
    it "orders longer overlaps before shorter overlaps" $
      let left  = OverlapCandidate (Fragment "AAA") (Fragment "BBB") 2
          right = OverlapCandidate (Fragment "CCC") (Fragment "DDD") 3
       in compareCandidates left right `shouldBe` LT

    it "prefers the lexicographically smaller prefix when lengths are equal" $
      let left  = OverlapCandidate (Fragment "ABC") (Fragment "XYZ") 3
          right = OverlapCandidate (Fragment "ABD") (Fragment "UVW") 3
       in compareCandidates left right `shouldBe` GT

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
