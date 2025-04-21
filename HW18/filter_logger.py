import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import logging

logging.basicConfig(filename='fastq_filter.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

class FastQFilter:
    def __init__(self, input_file: str):
        """
        Инициализация объектов FastQFilter

        Args:
            input_file (str): FASTQ-файл для фильтрации
        """
        self.input_file = input_file  
        
        
    def filter_by_length(self, min_length: int, max_length: int):
        """
        Фильтрует FASTQ-файл по длине

        Args:
            min_length (int): Минимальная длина
            max_length (int): Максимальная длина

        Returns:
            list: Список отфильтрованных записей
        """
        filtered_records = []
        try:
            for record in SeqIO.parse(self.input_file, "fastq"):
                if min_length <= len(record.seq) <= max_length:
                    filtered_records.append(record)
            logging.info(f"Фильтрация по длине завершена: {len(filtered_records)} записей удовлетворяют критериям.")
        except Exception as e:
            logging.error(f"Ошибка при фильтрации по длине: {e}")
        return filtered_records
    
    
    def filter_by_quality(self, min_quality: float):
        """
        Фильтрует FASTQ-файл по качеству
        Args:
            min_quality (float): Минимальное качество
        Returns:
            list: Список отфильтрованных записей
        """
        filtered_records = []
        try:
            for record in SeqIO.parse(self.input_file, "fastq"):
                quality_scores = record.letter_annotations["phred_quality"]
                avg_quality = sum(quality_scores) / len(quality_scores)
                if avg_quality >= min_quality:
                    filtered_records.append(record)
            logging.info(f"Фильтрация по качеству завершена: {len(filtered_records)} записей удовлетворяют критериям.")
        except Exception as e:
            logging.error(f"Ошибка при фильтрации по качеству: {e}")
        return filtered_records


    def filter_by_gc_content(self, min_gc: float, max_gc: float):
        """
        Фильтрует FASTQ-файл по GC-составу
        
        Args:
            min_gc (float): Минимальный GC-состав
            max_gc (float): Максимальный GC-состав
        Returns:
            list: Список отфильтрованных записей
        """
        filtered_records = []
        try:
            for record in SeqIO.parse(self.input_file, "fastq"):
                sequence = Seq(record.seq)
                gc_content = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
                if min_gc <= gc_content <= max_gc:
                    filtered_records.append(record)
            logging.info(f"Фильтрация по GC-содержанию завершена: {len(filtered_records)} записей удовлетворяют критериям.")
        except Exception as e:
            logging.error(f"Ошибка при фильтрации по GC-содержанию: {e}")
        return filtered_records
    
    
    def save_filtered_records(self, records, output_file: str):
        """
        Схраняет отфильтрованный FASTQ-файл
        
        Args:
            records (list): список отфильтрованных ридов
            output_file (str): путь до FASTQ-файла
     
        """
        try:
            with open(output_file, "w") as output_handle:
                SeqIO.write(records, output_handle, "fastq")
            logging.info(f'Отфильтрованные риды сохранены в файл: {output_file}')
        except Exception as e:
            logging.error(f'Ошибка при сохранении фильтрованных ридов {e}')               
        
try:
    parser = argparse.ArgumentParser(description = "Filter FASTQ-files by length, quality and GC-content")
    parser.add_argument("input_file", type = str, help = "Input FATSQ-file")
    parser.add_argument("output_file", type = str, help = "Output FASTQ-file")
    parser.add_argument("-mil", "--min_length", type = int, default = 0, help = "Minimum lenght of sequence")
    parser.add_argument("-mal", "--max_length", type = int, default = 1000, help = "Maximum length of sequence")
    parser.add_argument("-miq", "--min_quality", type = float, default = 0, help = "Minimum quality of sequence")
    parser.add_argument("-migc", "--min_gc", type = float, default = 0, help = "Minimum GC-content of sequence")
    parser.add_argument("-magc", "--max_gc", type = float, default = 100, help = "Maximum GC-content of sequence")
    args = parser.parse_args()
    
     # Создание экземпляра класса FastQFilter
    filter = FastQFilter(args.input_file)

        # Фильтрация по длине
    filtered_records = filter.filter_by_length(args.min_length, args.max_length)
    filtered_records = filter.filter_by_quality(args.min_quality)
    filtered_records = filter.filter_by_gc_content(args.min_gc, args.max_gc)
    filter.save_filtered_records(filtered_records, args.output_file)    

except SystemExit as e:
    print("Упс! Ошибка при запуске:", e)