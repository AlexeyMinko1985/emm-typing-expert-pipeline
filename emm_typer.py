import warnings
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from Bio import Align
import os

# --- БЛОК ПУТЕЙ (ОБЪЯВЛЕН В НАЧАЛЕ) ---
# Находим папку, в которой лежит скрипт
BASE = Path(__file__).resolve().parent

cons_dir = BASE / "results" / "consensuses_fastas"
db_dir = BASE / "data" / "db_cdc"
out_dir = BASE / "results" / "emm-types"
# --------------------------------------------------

"""
================================================================================
БИОИНФОРМАТИЧЕСКИЙ ПАЙПЛАЙН ЭКСПЕРТНОГО EMM-ТИПИРОВАНИЯ
================================================================================
ОПИСАНИЕ:
Данный скрипт автоматизирует процесс определения emm-типа Streptococcus pyogenes
на основе последовательностей гена emm. 

ЛОГИКА РАБОТЫ:
1. Загрузка баз данных CDC: trimmed (обрезанные до 180 б.п.) и untrimmed (полные).
2. Парсинг входных консенсусных последовательностей в формате FASTA.
3. Локальное выравнивание (Smith-Waterman) каждой пробы с референсами:
   - Match: +2, Mismatch: -3.
   - Определение % идентичности на основе длины перекрытия.
4. Валидация качества (QC):
   - Статус ПРИГОДЕН: Сходство >= 92% и длина чтения >= 180 б.п.
   - Статус НЕПРИГОДЕН: Несоответствие любому из критериев.
5. Генерация отчетности:
   - CSV-таблица со всеми параметрами выравнивания.
   - Визуальный отчет (Heatmap) со сводной статистикой прохождения контроля.
================================================================================
"""

warnings.filterwarnings("ignore", category=UserWarning)


def экспертное_типирование(query_seq, db_records):
    """Выполняет парное выравнивание и возвращает лучший хит."""
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1

    best_hit, max_score, best_idnt, best_len = "н/д", -1.0, 0.0, 0
    q_str = str(query_seq).upper()

    for hit_id, hit_seq in db_records.items():
        score = aligner.score(q_str, hit_seq.upper())
        if score > max_score:
            max_score = score
            best_hit = hit_id
            best_len = min(len(q_str), len(hit_seq))
            best_idnt = (score / (2 * best_len)) * 100

    return best_hit, round(min(best_idnt, 100.0), 2), best_len


def запустить_упорядоченный_пайплайн():
    # Автоматическая генерация даты для отчета
    today = datetime.now().strftime('%d.%m.%Y')

    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"--- ЗАПУСК ТИПИРОВАНИЯ (ДАТА: {today}) ---")

    try:
        db_t = {rec.id: str(rec.seq) for rec in SeqIO.parse(db_dir / "alltrimmed.tfa", "fasta")}
        db_u = {rec.id: str(rec.seq) for rec in SeqIO.parse(db_dir / "alluntrimmed.tfa", "fasta")}
    except Exception as e:
        print(f"❌ ОШИБКА ЗАГРУЗКИ БАЗ: {e}")
        return

    raw_files = []
    for f in cons_dir.glob('*'):
        if f.is_file() and f.name.lower().startswith('consensus'):
            match = re.search(r'(\d+)', f.name)
            if match: raw_files.append((int(match.group(1)), f))

    raw_files.sort(key=lambda x: x[0])

    report_data = []
    for s_id, f_path in raw_files:
        try:
            rec = SeqIO.read(f_path, "fasta")
            t_id, t_idnt, t_len = экспертное_типирование(rec.seq, db_t)
            u_id, _, _ = экспертное_типирование(rec.seq, db_u)

            def clean_emm(name):
                clean_name = re.sub(r'emm|em|\.0', '', name, flags=re.IGNORECASE)
                return f"emm{clean_name}"

            verdict_full = f"{clean_emm(t_id)}.0({clean_emm(u_id)}.0)"

            if t_idnt >= 92.0 and t_len >= 180:
                qc_val = 1.0
                status_text = "ПРИГОДЕН"
            else:
                qc_val = 0.0
                status_text = "НЕПРИГОДЕН"

            report_data.append({
                'ID': s_id, 'Результат': verdict_full,
                'Идентификация_%': t_idnt, 'Перекрытие_bp': t_len,
                'Статус_QC': qc_val, 'Вердикт': status_text
            })
            print(f"✅ Обработан ID {s_id}: {verdict_full}")
        except Exception as e:
            print(f"❌ Ошибка в файле {f_path.name}: {e}")

    if report_data:
        df = pd.DataFrame(report_data)
        total = len(df)
        passed = len(df[df['Статус_QC'] == 1.0])
        failed_list = df[df['Статус_QC'] == 0.0]['ID'].tolist()

        # Разделение результата на Тип и Подтип для визуализации
        df['Тип'] = df['Результат'].apply(lambda x: x.split('(')[0] if '(' in x else x)
        df['Подтип'] = df['Результат'].apply(lambda x: x.split('(')[1].strip(')') if '(' in x else "н/д")

        # Набор колонок для тепловой карты (Твои обозначения сохранены)
        display_cols = ['Тип', 'Подтип', 'Идентификация_%', 'Перекрытие_bp', 'Вердикт']
        plot_df_numeric = df.set_index('ID')[
            ['Статус_QC', 'Статус_QC', 'Идентификация_%', 'Перекрытие_bp', 'Статус_QC']]
        annot_matrix = df.set_index('ID')[display_cols].astype(str)

        # Визуализация
        plt.figure(figsize=(14, 6 + (len(df) * 0.3)))
        ax = sns.heatmap(plot_df_numeric, annot=annot_matrix, fmt="", cmap="RdYlGn", center=0.5, cbar=False)

        # ПЕРЕНОС ЗАГОЛОВКОВ НАВЕРХ (Без изменения названий)
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.set_xticklabels(['Тип', 'Подтип', 'Идентичность_%', 'Перекрытие_bp', 'Вердикт'], fontsize=10)

        plt.title(f"ОТЧЕТ ПО ЭКСПЕРТНОМУ ТИПИРОВАНИЮ ({today})", fontsize=16, fontweight='bold', pad=40)

        # Инфо-панель (Сохранена полностью)
        stats_info = (
            f"ОБЩАЯ СТАТИСТИКА:\n"
            f"━━━━━━━━━━━━━━━━━━━━\n"
            f"Всего проанализировано: {total}\n"
            f"✅ Прошли контроль: {passed}\n"
            f"❌ Не прошли контроль: {total - passed}\n"
            f"━━━━━━━━━━━━━━━━━━━━\n"
            f"ID НЕПРИГОДНЫХ:\n"
            f"{', '.join(map(str, failed_list)) if failed_list else 'отсутствуют'}"
        )

        plt.figtext(0.80, 0.5, stats_info, fontsize=10, family='monospace', va='center',
                    bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round'))

        plt.subplots_adjust(left=0.1, right=0.75, bottom=0.1, top=0.85)

        save_png = out_dir / f"Expert_emm_Typing_Report_{today}.png"
        plt.savefig(save_png, dpi=300, bbox_inches='tight')
        df.to_csv(out_dir / f"Final_emm_Typing_Results_{today}.csv", index=False, encoding='utf-8-sig')

        print(f"\n📊 АНАЛИЗ ЗАВЕРШЕН. Результаты в: {out_dir}")
    else:
        print("❌ Ошибка: Нет данных для анализа.")


if __name__ == "__main__":
    # Переходим в папку скрипта, чтобы относительные пути в начале сработали
    os.chdir(BASE)
    запустить_упорядоченный_пайплайн()

