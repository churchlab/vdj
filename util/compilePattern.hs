-- This code takes bitmask combs and generates corresponding
-- C code for applying the comb

import Text.Printf
import Data.List 
import System.Environment

-- count, mask, id -> count, mask, shift
computeShift :: [(Int, Int, Bool)] -> [(Int, Int, Int)]
computeShift al =
  snd $ foldr shift (0, []) al
  where shift (c, m, i) (o, l) = 
          if i
          then (o, ((c, m, o):l))
          else (o + c, l)

-- count, offset, shift -> mask, shift
makeShiftMask :: [(Int, Int, Int)] -> [(Int,Int)]
makeShiftMask al =
  map (\ (c, offset, shift) -> ((foldl1 bin2dec [1 | x <- [1..c]]) * (2 ^ offset), shift)) al 
  where bin2dec a b = a * 2 + b

count :: (Eq a) => [a] -> [(Int, a)]
count al = 
  foldr squish [] $ map (\ a -> (1,a)) al

squish :: (Eq a) => (Int, a) -> [(Int, a)] -> [(Int, a)]
squish a [] = [a]
squish (x,a) l@((y,b):tl) = 
  if a == b
  then ((x + y), a):tl
  else (x,a):l

 -- count, id -> count, mask, id
makeOffset :: [(Int, Bool)] -> [(Int, Int, Bool)]
makeOffset a = 
  snd $ foldr (\ (i, a) (n, r) -> (n + i, (i, n, a):r)) (0, []) a

maskShiftToC :: String -> (Int, Int) -> String
maskShiftToC n (mask, shift) =
  printf "((0x%x & %s) >> %d)" mask n shift

binToNum :: [Bool] -> String
binToNum n = binToNumH n 0
binToNumH :: [Bool] -> Int -> String
binToNumH [] n = printf "0x%X" n
binToNumH (True:tl)  n = binToNumH tl (2 * n + 1)
binToNumH (False:tl) n = binToNumH tl (2 * n)

patternToMasks :: (Show a) => a -> [Bool] -> String
patternToMasks x p =
 let sm = makeShiftMask $ computeShift $ makeOffset $ count p in
   let accShift = concat $ intersperse " | " $ map (maskShiftToC ("xx" ++ (show x))) sm
       mskShift = "(" ++ (binToNum p) ++ " & mask)" in
     concat $ "if(!(" : mskShift : ") && (ii >  " : (show (length p)) : ")) { insertBump(hm, " : accShift : "); }" : []

processComb :: (Show a) => a -> String -> String
processComb var c =
  patternToMasks var $ map (\ x -> if x == '1' then True else False) c 

main = do
  n <- getArgs
  putStr $ unlines $ zipWith ($) (map processComb "xx") n
  
