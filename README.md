# Программный комплекс для автоматизации реконструкции, сборки консенсусных последовательностей и emm-типирования Streptococcus pyogenes

Данный проект представляет собой локальную автоматизированную реализацию протоколов CDC по молекулярному субтипированию стрептококков группы А. Логика пайплайна полностью соответствует техническим стандартам "Streptococci Group A Subtyping Request Form" и алгоритмам "Blast 2.0 Server".

## Научная и методологическая база
Разработка базируется на отраслевых стандартах и верифицированных базах данных следующих подразделений:
*   U.S. Centers for Disease Control and Prevention (CDC)
*   Division of Bacterial Diseases
*   Biotechnology Core Facility Computing Laboratory

## Описание системы
Комплекс автоматизирует два ключевых этапа биоинформатического анализа:
1. Реконструкция, сборка и верификация консенсусных последовательностей из первичных данных хроматограмм (Sanger sequencing).
2. Экспертное emm-типирование на основе референсных баз данных CDC.

## Функциональные возможности

### 1. Модуль реконструкции, сборки и верификации (Sanger Assembly QC)
*   Автоматическая идентификация пар ридов (Forward/Reverse) по цифровому идентификатору образца.
*   Качественный тримминг: автоматическое отсечение низкокачественных фрагментов последовательности (порог Phred Score > 20).
*   Алгоритм N-коррекции: интеллектуальное устранение неопределенных нуклеотидов в зонах перекрытия ридов на основе сравнительного анализа Phred-баллов.
*   Глобальное выравнивание (PairwiseAligner) для обеспечения достоверности верифицированного консенсуса.

### 2. Модуль экспертного типирования (emm-typing)
*   Локальное выравнивание по алгоритму Смита-Ватермана (полная совместимость с Blast 2.0 Server).
*   Двойная верификация по референсным базам CDC (trimmed и untrimmed).
*   Валидация результатов по критериям Biotechnology Core Facility (CDC): идентичность не менее 92%, длина перекрытия не менее 150 п.н.
*   Автоматическое разделение результатов на Тип и Подтип.

## Модуль визуализации данных (Seaborn Data Visualization)
Для графического представления результатов анализа в систему интегрирован модуль визуализации на базе библиотеки Seaborn. Это позволяет формировать отчеты экспертного уровня:
*   **Автоматическая генерация тепловых карт (Heatmaps):** Наглядное представление результатов типирования с разделением на типы и подтипы.
*   **Цветовое кодирование (Color Mapping):** Применение палитр для мгновенной идентификации качества данных и вердикта QC.
*   **Динамическая статистика:** Вывод ключевых метрик (количество групп, собранных пар, статистика брака) непосредственно на графические отчеты.

### Сгенерированные отчеты:
1. **Отчет о качестве сборки консенсусов:** 
![QC Report](Consensus_Quality_Report.png)

2. **Отчет об экспертном emm-типировании:** 
![Typing Report](Expert_emm_Typing_Report.png)

## Структура проекта (Автоматическая настройка)
Пайплайн спроектирован для работы с минимальной предварительной настройкой. Все необходимые рабочие директории создаются скриптами автоматически при первом запуске.

*   **data/db_cdc/** — Директория для референсных баз данных CDC. Заполняется пользователем.
*   **data/raw_input/** — (Автоматическое создание) Директория для размещения входящих файлов хроматограмм (.ab1).
*   **results/consensuses_fastas/** — (Автоматическое создание) Директория для сохранения верифицированных консенсусов.
*   **results/emm-types/** — (Автоматическое создание) Директория для итоговых экспертных отчетов (PNG и CSV).
*   **requirements.txt** — Список необходимых библиотек для развертывания программной среды.

## Инструкция по эксплуатации
1. Установка необходимых библиотек: `pip install -r requirements.txt`.
2. Размещение исходных файлов .ab1 в автоматически созданной директории `data/raw_input/`.
3. Последовательный запуск модулей сборки и типирования.

---
Разработчик: Алексей Минко (Microbiologist & Bioinformatics Developer).
Методология: Biotechnology Core Facility Computing Laboratory (CDC).




---

# Automated Software Package for Sanger Consensus Assembly and emm-Typing of Streptococcus pyogenes

This project is a localized automated implementation of CDC protocols for molecular subtyping of Group A Streptococci. The pipeline logic fully complies with the "Streptococci Group A Subtyping Request Form" standards and "Blast 2.0 Server" algorithms.

## Scientific and Methodological Background
The development is based on industry standards and verified databases from the following institutions:
*   U.S. Centers for Disease Control and Prevention (CDC)
*   Division of Bacterial Diseases
*   Biotechnology Core Facility Computing Laboratory

## System Overview
The package automates two key stages of bioinformatic analysis:
1. Reconstruction, assembly, and verification of consensus sequences from raw Sanger sequencing chromatograms.
2. Expert emm-typing based on reference CDC databases.

## Functional Features

### 1. Reconstruction, Assembly, and Verification Module (Sanger Assembly QC)
*   Automatic identification of read pairs (Forward/Reverse) by sample ID.
*   Quality Trimming: Automatic removal of low-quality sequence fragments (Phred Score threshold > 20).
*   N-Correction Algorithm: Intelligent elimination of uncertain nucleotides in overlap zones based on Phred score comparison.
*   Global Alignment (PairwiseAligner) to ensure the reliability of the verified consensus.

### 2. Expert emm-Typing Module
*   Local alignment using the Smith-Waterman algorithm (full compatibility with Blast 2.0 Server).
*   Double verification using CDC reference databases (trimmed and untrimmed).
*   Result validation according to Biotechnology Core Facility (CDC) criteria: identity ≥ 92%, overlap length ≥ 150 bp.
*   Automatic separation of results into Type and Subtype.

## Data Visualization Module (Seaborn)
The system integrates a visualization module based on the Seaborn library for expert-level reporting:
*   **Automated Heatmaps:** Visual representation of typing results with clear type/subtype separation.
*   **Color Mapping:** Palette application for instant identification of data quality and QC verdicts.
*   **Dynamic Statistics:** Key metrics (group count, assembled pairs, failure statistics) displayed directly on graphical reports.

## Project Structure (Automated Setup)
The pipeline is designed for "Zero-Config" operation. All necessary working directories are created automatically by scripts upon the first run.

*   **data/db_cdc/** — Directory for CDC reference databases. (User-provided).
*   **data/raw_input/** — (Auto-created) Directory for incoming chromatogram files (.ab1).
*   **results/consensuses_fastas/** — (Auto-created) Directory for verified consensus storage.
*   **results/emm-types/** — (Auto-created) Directory for final expert reports (PNG and CSV).
*   **requirements.txt** — List of dependencies for environment deployment.

## Installation and Usage
1. Install dependencies: `pip install -r requirements.txt`.
2. Place raw .ab1 files into the automatically created `data/raw_input/` directory.
3. Run assembly and typing modules sequentially.

---
**Developer:** Alexey Minko (Microbiologist & Bioinformatics Developer).
**Methodology:** Biotechnology Core Facility Computing Laboratory (CDC).

