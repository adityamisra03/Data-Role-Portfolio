-- Exploratory Data Analysis

-- CTE
-- Rolling Total
-- Window Function

-- GOAL: Rank the top 5 companies that laid off the most people each year

SELECT *
FROM layoffs_staging2;

SELECT MAX(total_laid_off), MAX(percentage_laid_off)
FROM layoffs_staging2;
-- Max number of people laid off, (12000 people and  100% people laid off)

SELECT *
FROM layoffs_staging2
WHERE percentage_laid_off = 1
ORDER BY total_laid_off DESC;
-- highest 100% laid off by number

SELECT *
FROM layoffs_staging2
WHERE percentage_laid_off = 1
ORDER BY funds_raised_millions DESC;
-- highest 100% laid off by millions raised

SELECT company, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company
ORDER BY 2 DESC;
-- Group by company ordered by SUM(total_laid_off)

SELECT MIN(`date`), MAX(`date`)
FROM layoffs_staging2;
-- almost 3 years of data for layoffs

SELECT industry, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY industry
ORDER BY 2 DESC;
-- Group by industry ordered by SUM(total_laid_off)

SELECT YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY YEAR(`date`)
ORDER BY 1 DESC;
-- Around 125,700 layoffs in 2023 and 160661 layoffs in 2022

SELECT *
FROM layoffs_staging2;

SELECT stage, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY stage
ORDER BY 1 DESC;

SELECT stage, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY stage
ORDER BY 2 DESC;

SELECT company, AVG(percentage_laid_off)
FROM layoffs_staging2
GROUP BY company
ORDER BY 2 DESC;
-- Not that useful

SELECT SUBSTRING(`date`,6,2) AS `MONTH`, SUM(total_laid_off) -- Month
FROM layoffs_staging2
GROUP BY `MONTH`
ORDER BY 1 ASC;

SELECT SUBSTRING(`date`,1,7) AS `MONTH`, SUM(total_laid_off) -- Year and Month
FROM layoffs_staging2
WHERE SUBSTRING(`date`,1,7) IS NOT NULL
GROUP BY `MONTH`
ORDER BY 1 ASC;

SELECT *
FROM layoffs_staging2;

WITH Rolling_Total AS -- CTE
(
SELECT SUBSTRING(`date`,1,7) AS `MONTH`, SUM(total_laid_off) AS total_off -- Year and Month
FROM layoffs_staging2
WHERE SUBSTRING(`date`,1,7) IS NOT NULL
GROUP BY `MONTH`
ORDER BY 1 ASC
)
SELECT `MONTH`, total_off,
SUM(total_off) OVER(ORDER BY `MONTH`) AS rolling_total -- Window Function
FROM Rolling_Total;
-- See the progression of layoffs through rolling total

SELECT company, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company
ORDER BY 2 DESC;

SELECT company, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company, YEAR(`date`)
ORDER BY company ASC;
-- Ordered by company name and year laidoff

SELECT company, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company, YEAR(`date`)
ORDER BY 3 DESC;
-- Want to see using CTE

WITH Company_Year (company, years, total_laid_off) AS -- CTE
(
SELECT company, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company, YEAR(`date`)
)
SELECT *, DENSE_RANK() OVER (PARTITION BY years ORDER BY total_laid_off DESC)
FROM Company_Year
WHERE years IS NOT NULL; -- some null data messing up ranking
-- Will go through data and rank them without repeats as separate (DENSE_RANK)
-- order by total_laid_off split up and ranking from each separate year

SELECT company, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company, YEAR(`date`)
ORDER BY 3 DESC;
-- Want to see using CTE

WITH Company_Year (company, years, total_laid_off) AS -- CTE
(
SELECT company, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company, YEAR(`date`)
)
SELECT *, DENSE_RANK() OVER (PARTITION BY years ORDER BY total_laid_off DESC) AS Ranking
FROM Company_Year
WHERE years IS NOT NULL -- some null data messing up ranking
ORDER BY Ranking ASC;
-- Will go through data and rank them without repeats as separate (DENSE_RANK)
-- order by total_laid_off grouping all years together on the ranking (1's all together, 2's all together, etc.)




WITH Company_Year (company, years, total_laid_off) AS -- CTE and subqueries
(
SELECT company, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company, YEAR(`date`)
), Company_Year_Rank AS
(
SELECT *, DENSE_RANK() OVER (PARTITION BY years ORDER BY total_laid_off DESC) AS Ranking
FROM Company_Year
WHERE years IS NOT NULL -- some null data messing up ranking
)
SELECT *
FROM Company_Year_Rank
WHERE Ranking <= 5;
-- ranking total_laid_off of the top 5 companies each year that laid off people each year
-- 1: Uber 2020, Bytedance 2021, Meta 2022, Google 2023