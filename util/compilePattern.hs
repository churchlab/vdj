-- This code takes bitmask combs and generates corresponding
-- C code for applying the comb

import Text.Printf
import Data.List 
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

maskShiftToC :: (Int, Int) -> String
maskShiftToC (mask, shift) =
  printf "((0x%x & xx) >> %d)" mask shift

patternToMasks :: [Bool] -> String
patternToMasks p =
  (concat $ "return " : (intersperse " | " $ map maskShiftToC $ makeShiftMask $ computeShift $ makeOffset $ count p)) ++ ";"

processComb :: String -> String
processComb c =
  patternToMasks $ map (\ x -> if x == '1' then True else False) c 
