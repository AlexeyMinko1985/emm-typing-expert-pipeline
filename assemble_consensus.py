import os
import re
import warnings
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from datetime import datetime
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
# Принудительно переходим в папку, где лежит сам скрипт
os.chdir(Path(__file__).resolve().parent)

# Теперь, когда мы уже в нужной папке, пути пишем просто:
PROJECT_ROOT = Path(".").resolve()
DATA_DIR = PROJECT_ROOT / "data"
RAW_INPUT_DIR = DATA_DIR / "raw_input"
RESULTS_DIR = PROJECT_ROOT / "results"
CONSENSUS_DIR = RESULTS_DIR / "consensuses_fastas"

# ----------------------------------------------

"""
================================================================================
БИОИНФОРМАТИЧЕСКИЙ ПАЙПЛАЙН СБОРКИ КОНСЕНСУСОВ (ОТЧЕТ ОТ 20.01.2026)
================================================================================
ФУНКЦИОНАЛ:
1. ВСЕЯДНОСТЬ: Читает .ab1, .fasta, .txt, .seq.
2. ПОИСК ПАР: Сначала ищет цифру (ID), затем метки направления (F, R, FWD, REV).
3. АВТО-РЕВЕРС: Делает Reverse Complement для обратных последовательностей.
4. КАЧЕСТВО: Сравнивает Phred баллы и собирает консенсус по лучшим значениям.
5. ЛЕГЕНДА: Графические цветовые полоски, статистика и вердикт в правом блоке.
================================================================================
"""

warnings.filterwarnings("ignore", category=UserWarning)


def получить_статус(q_score):
    """Определяет категорию качества на основе среднего Phred балла."""
    if q_score >= 30:
        return "Отличный"
    elif q_score >= 20:
        return "Посредственный"
    else:
        return "Непригодный"


def тримминг_качества(seq, qual, threshold=20):
    """НОВЫЙ БЛОК: Отсекает концы рида, где Phred баллы ниже порога (threshold)."""
    start, end = 0, len(qual)
    for i, q in enumerate(qual):
        if q >= threshold:
            start = i
            break
    for i, q in enumerate(reversed(qual)):
        if q >= threshold:
            end = len(qual) - i
            break
    return seq[start:end], qual[start:end]


def загрузить_данные(file_path):
    """Загружает последовательность и качество из файла."""
    ext = file_path.suffix.lower()
    try:
        if ext in ['.ab1', '.abi', '.abif']:
            trace = SeqIO.read(file_path, "abi")
            s, q = str(trace.seq).upper(), list(trace.letter_annotations.get("phred_quality"))
            return тримминг_качества(s, q)  # ПРИМЕНЕН ТРИММИНГ ДЛЯ AB1
        else:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            seq_raw = "".join(content.splitlines()[1:]) if content.startswith('>') else content
            clean_seq = "".join(re.findall(r'[atgcnATGCN]', seq_raw)).upper()
            return clean_seq, [40] * len(clean_seq)  # Исправлен баг с *
    except Exception:
        return None, None


def запустить_пайплайн(input_path):
    today = datetime.now().strftime("%d.%m.%Y")
    src_dir = Path(input_path).resolve()
    
    # СОЗДАЕМ ТОЛЬКО RESULTS
    CONSENSUS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"--- ЗАПУСК СБОРКИ ({today}) ---")

    valid_exts = ('.ab1', '.abi', '.abif', '.fasta', '.fa', '.seq', '.txt')
    all_files = [f for f in src_dir.rglob('*') if f.suffix.lower() in valid_exts]

    # Настройка выравнивателя
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1

    # 1. ГРУППИРОВКА ПО ЦИФРАМ (ID)
    strain_groups = {}
    for f_path in all_files:
        match = re.search(r'\d+', f_path.name)
        if match:
            s_id = match.group()
            if s_id not in strain_groups: strain_groups[s_id] = []
            strain_groups[s_id].append(f_path)

    report_data, unpaired_files = [], []
    FWD_L = ['FOR', 'FWD', 'DIR', '_F', 'F.', 'F_']
    REV_L = ['REV', 'RV', 'SR', '_R', 'R.', 'R_']

    # 2. ПОИСК НАПРАВЛЕНИЙ ВНУТРИ ГРУППЫ
    for s_id in sorted(strain_groups.keys(), key=lambda x: int(x)):
        files = strain_groups[s_id]
        f_match = [f for f in files if any(l in f.name.upper() for l in FWD_L) or f.stem.upper().endswith('F')]
        r_match = [f for f in files if any(l in f.name.upper() for l in REV_L) or f.stem.upper().endswith('R')]

        if f_match and r_match:
            fp, rp = f_match[0], r_match[0]
        elif len(files) == 2:
            fp, rp = files[0], files[1]
        else:
            for f in files: unpaired_files.append(f.name)
            continue

        try:
            s1, q1 = загрузить_данные(fp)
            s2_raw, q2_raw = загрузить_данные(rp)
            if s1 and s2_raw:
                s2 = str(Seq(s2_raw).reverse_complement())
                q2 = q2_raw[::-1]

                # Внедрение выравнивания
                alignments = aligner.align(s1, s2)
                best_aln = alignments[0]
                aln_s1, aln_s2 = best_aln[0], best_aln[1]

                consensus, quals = [], []
                idx1, idx2 = 0, 0

                for char1, char2 in zip(aln_s1, aln_s2):
                    v1 = q1[idx1] if char1 != '-' else 0
                    v2 = q2[idx2] if char2 != '-' else 0

                    if char1 == char2 and char1 != '-':
                        consensus.append(char1); quals.append(max(v1, v2))
                    elif char1 == '-':
                        consensus.append(char2); quals.append(v2)
                    elif char2 == '-':
                        consensus.append(char1); quals.append(v1)
                    elif char1 == 'N' and char2 not in ['N', '-']:
                        consensus.append(char2); quals.append(v2)
                    elif char2 == 'N' and char1 not in ['N', '-']:
                        consensus.append(char1); quals.append(v1)
                    else:
                        if v1 >= v2:
                            consensus.append(char1); quals.append(v1)
                        else:
                            consensus.append(char2); quals.append(v2)

                    if char1 != '-': idx1 += 1
                    if char2 != '-': idx2 += 1

                avg_q = sum(quals) / len(quals) if quals else 0
                report_data.append({'ID': s_id, 'Phred_Q': avg_q, 'Статус': получить_статус(avg_q)})

                rec = SeqRecord(Seq("".join(consensus)), id=f"ID_{s_id}", description=f"Q:{avg_q:.1f}")
                SeqIO.write(rec, CONSENSUS_DIR / f"Consensus_{s_id}.fasta", "fasta")
                print(f"✅ ID {s_id}: Собрано (Тримминг и N-коррекция применены)")
        except Exception as e:
            print(f"❌ Ошибка в ID {s_id}: {e}")

    # 3. ВИЗУАЛИЗАЦИЯ
    if report_data:
        df = pd.DataFrame(report_data)
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(figsize=(max(12, len(df) * 0.6), 8))

        palette = {"Отличный": "#2ecc71", "Посредственный": "#f1c40f", "Непригодный": "#e74c3c"}
        sns.barplot(data=df, x='ID', y='Phred_Q', hue='Статус', palette=palette, ax=ax, dodge=False)

        ax.set_ylim(0, 50)
        ax.set_ylabel('Качество (Phred баллы Q)', fontsize=14, fontweight='bold', labelpad=25)
        ax.set_xlabel('Идентификатор штамма (ID)', fontsize=12, fontweight='bold')
        ax.set_title(f'АНАЛИЗ КАЧЕСТВА КОНСЕНСУСОВ — {today}', fontsize=16, pad=30, fontweight='bold')

        counts = df['Статус'].value_counts()
        legend_elements = [
            Line2D([0], [0], color='#2ecc71', lw=8, label=f'Отличный (Q30+): {counts.get("Отличный", 0)}'),
            Line2D([0], [0], color='#f1c40f', lw=8, label=f'Посредственный (Q20-30): {counts.get("Посредственный", 0)}'),
            Line2D([0], [0], color='#e74c3c', lw=8, label=f'Непригодный (Q<20): {counts.get("Непригодный", 0)}')
        ]
        ax.legend(handles=legend_elements, title="ЛЕГЕНДА", bbox_to_anchor=(1.02, 1), loc='upper left')

        stats_text = (f"ОБЩАЯ СТАТИСТИКА:\n{'-' * 20}\nВсего групп: {len(strain_groups)}\n"
                      f"Собрано пар: {len(df)}\nБез пары: {len(unpaired_files)}")

        plt.subplots_adjust(left=0.15, right=0.75, bottom=0.15, top=0.85)
        fig.text(0.76, 0.4, stats_text, fontsize=10, family='monospace', va='top',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='#dee2e6', boxstyle='round'))

        # ПУТЬ СОХРАНЕНИЯ ГРАФИКА
        save_path = RESULTS_DIR / f"Consensus_Quality_Report_{today}.png"
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"--- ГОТОВО. ОТЧЕТ: {save_path} ---")
        plt.show()
    else:
        print("❌ Ошибка: Данные не найдены.")


if __name__ == "__main__":
    if not RAW_INPUT_DIR.exists():
        print(f"❌ Ошибка: Папка не найдена: {RAW_INPUT_DIR}")
    else:
        запустить_пайплайн(RAW_INPUT_DIR)
