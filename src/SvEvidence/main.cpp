#include "ToolBase.h"
#include "BamWriter.h"
#include "Helper.h"
#include <QFile>
#include <QTextStream>
#include "BasicStatistics.h"
#include <QElapsedTimer>
#include <QDebug>

class ConcreteTool
		: public ToolBase
{
	Q_OBJECT

public:
	ConcreteTool(int& argc, char *argv[])
		: ToolBase(argc, argv)
	{
	}

	virtual void setup()
	{
		setDescription("Computes some statistics of the reads in different areas of the genome.");
		addInfile("in", "Input BAM file.", false);
		addOutfile("out", "Output CSV file containing the computed statistics.", false);

		// optional
		addInt("binSize", "Defines the bin size in which the genome is divided.", true, 100);
		addInt("qFilter", "Determines the minimum quality score for each alignment.", true, 0);
		addInt("ISDistSize", "Determines the number of alignments used for the insert size distribution.", true, 100000);

		// changelog
		changeLog(2019, 4, 17, "Initial version of this tool");
		changeLog(2019, 4, 18, "Added insert size distribution");
		changeLog(2019, 4, 23, "Improved insert size distribution statistics");
	}

	void compute_insert_size_cutoffs(){
		//compute insert size distribution
		BamReader reader(getInfile("in"));
		QVector<double> insert_size_distribution = QVector<double>(insertSizeDistributionSize);
		//process alignments
		int alignment_counter = 0;
		BamAlignment al;
		while (reader.getNextAlignment(al)) {

			// skip not usefull reads
			if (al.isSecondaryAlignment()) continue;
			if (al.isUnmapped()) continue;
			if (al.mappingQuality() < minimumQualityScore) continue;


			// store insert size of this alignment
			insert_size_distribution[alignment_counter] = std::abs(al.insertSize());

			alignment_counter++;
			if(alignment_counter >= insertSizeDistributionSize){
				break;
			}
		}

		// exeption if BAM file is to small
		if (alignment_counter < insertSizeDistributionSize){
			// bam file contains fewer elements than distribution vector:
			//TODO: Error handling

			// remove not used elements of the qvector
			for (int index=alignment_counter; index < insertSizeDistributionSize; index++){
				insert_size_distribution.removeLast();
			}
		}

		// compute robust statistics
		// sort qVector
		std::sort(insert_size_distribution.begin(), insert_size_distribution.end());
		double median = BasicStatistics::median(insert_size_distribution);
		double mad = BasicStatistics::mad(insert_size_distribution, median);

		//filter by median and mad:
		double lower_cutoff = median - 3 * mad;
		double upper_cutoff = median + 3 * mad;
		QVector<double> filtered_insert_size_distribution = QVector<double>();
		foreach (double value, insert_size_distribution) {
			if((value >= lower_cutoff) && (value <= upper_cutoff)){
				filtered_insert_size_distribution.append(value);
			}
		}

		// compute mean and stddev on filtered data
		double mean = BasicStatistics::mean(filtered_insert_size_distribution);
		double stddev = BasicStatistics::stdev(filtered_insert_size_distribution, mean);

		// store global cutoffs
		minInsertSize = mean - 3 * stddev;
		maxInsertSize = mean + 3 * stddev;

	}

	QString get_time_string(qint64 milliseconds){
		QTime time(0,0,0);
		time = time.addMSecs(milliseconds);
		return time.toString("hh:mm:ss.zzz");
	}




	virtual void main()
	{
		//init
		QTextStream out(stdout);

		binSize = getInt("binSize");
		insertSizeDistributionSize = getInt("ISDistSize");
		minimumQualityScore = getInt("qFilter");

		BamReader reader(getInfile("in"));

		QElapsedTimer timer;
		timer.start();


		//compute insert size distribution
		out << "calculating insert size cutoff values based on the first "
			<< insertSizeDistributionSize << " reads ..." << endl;

		compute_insert_size_cutoffs();
		// report cutoff
		out << "insert size distribution: " << endl;
		out << "lower cutoff: " << minInsertSize << endl;
		out << "upper cutoff: " << maxInsertSize << endl;
		qint64 insert_size_cutoff_calculation_runtime = timer.restart();
		out << "(calculation took " << get_time_string(insert_size_cutoff_calculation_runtime)
			<< ")" << endl << endl;


		out << "parsing BAM file ..." << endl;
		// setup bins
		QList<Chromosome> chromosomes = reader.chromosomes();
		QVector<QVector<int>> read_counts =
				QVector<QVector<int>>(chromosomes.size()); //contains the read counts for each bin
		QVector<QVector<long long>> total_bases =
				QVector<QVector<long long>>(chromosomes.size()); //contains the total number of bases mapped to this bin
		QVector<QVector<long long>> soft_clipped_bases =
				QVector<QVector<long long>>(chromosomes.size()); //contains the number of soft-clipped bases for each bin
		QVector<QVector<int>> lower_insert_size_outliers =
				QVector<QVector<int>>(chromosomes.size()); //contains the number of reads which have a too small insert size
		QVector<QVector<int>> upper_insert_size_outliers =
				QVector<QVector<int>>(chromosomes.size()); //contains the number of reads which have a too large insert size
		QVector<QVector<int>> forward_orientated_reads =
				QVector<QVector<int>>(chromosomes.size()); //contains the number of ">>" orientated reads for each bin
		QVector<QVector<int>> reverse_orientated_reads =
				QVector<QVector<int>>(chromosomes.size()); //contains the number of "<<" orientated reads for each bin
		QVector<QVector<int>> correct_orientated_reads =
				QVector<QVector<int>>(chromosomes.size()); //contains the number of "><" orientated reads for each bin
		QVector<QVector<int>> switched_orientated_reads =
				QVector<QVector<int>>(chromosomes.size()); //contains the number of "<>" orientated reads for each bin


		// initialize bins
		for(int i=0; i < chromosomes.size(); i++){
			int n_bins = reader.chromosomeSize(chromosomes[i]) / binSize + 1; //compute the number of bins for each chromosome
			out << reader.chromosome(i).str() << ": " << n_bins << endl;
			read_counts[i] =  QVector<int>(n_bins, 0);
			total_bases[i] = QVector<long long>(n_bins, 0);
			soft_clipped_bases[i] = QVector<long long>(n_bins, 0);
			upper_insert_size_outliers[i] = QVector<int>(n_bins, 0);
			lower_insert_size_outliers[i] = QVector<int>(n_bins, 0);
			forward_orientated_reads[i] = QVector<int>(n_bins, 0);
			reverse_orientated_reads[i] = QVector<int>(n_bins, 0);
			correct_orientated_reads[i] = QVector<int>(n_bins, 0);
			switched_orientated_reads[i] = QVector<int>(n_bins, 0);
		}



		int parsed_alignments = 0;
		int skipped_alignments = 0;
		BamAlignment al;
		while (reader.getNextAlignment(al))
		{
			skipped_alignments++;

			// skip not usefull reads
			if (al.isSecondaryAlignment()) continue;
			if (al.isUnmapped()) continue;
			if (al.mappingQuality() < minimumQualityScore) continue;

			skipped_alignments--;
			parsed_alignments++;

			const int chr_id = al.chromosomeID();
			const int start_pos = al.start();
			//const int end_pos = al.end();

			//determine bin:
			int bin = start_pos / binSize;


			// count overall reads
			read_counts[chr_id][bin]++;

			// count soft-clipped bases
			const QList<CigarOp> cigar_data = al.cigarData();
			foreach(const CigarOp& op, cigar_data)
			{
				if (op.Type==BAM_CSOFT_CLIP)
				{
					soft_clipped_bases[chr_id][bin] += op.Length;

				}else if (op.Type==BAM_CMATCH ||
						  op.Type==BAM_CEQUAL ||
						  op.Type==BAM_CDIFF) {
					// count overall aligned bases
					total_bases[chr_id][bin] += op.Length;

				}
			}

			if (al.isPaired()){
				// check for insert size outliers
				int abs_insertSize = std::abs(al.insertSize());
				if(abs_insertSize > maxInsertSize){
					upper_insert_size_outliers[chr_id][bin]++;
				}
				if(abs_insertSize < minInsertSize){
					lower_insert_size_outliers[chr_id][bin]++;
				}

				//check orientation:
				if (al.isReverseStrand() == al.isMateReverseStrand()){
					// both reads have the same orientation
					if (al.isMateReverseStrand()){
						//both reads are reversed
						reverse_orientated_reads[chr_id][bin]++;
					}else{
						// both reads are forward orientated
						forward_orientated_reads[chr_id][bin]++;
					}
				}else{
					// reads have different orientations
					if (al.isReverseStrand() == (al.start() > al.mateStart())){
						// reads have correct orientation
						correct_orientated_reads[chr_id][bin]++;
					}else{
						// reads have switched orientation
						switched_orientated_reads[chr_id][bin]++;
					}
				}
			}

			if(parsed_alignments%1000000 == 0){
				out << parsed_alignments / 1000000 << "M alignments parsed" << endl;
			}
		}
		out << parsed_alignments << " alignments parsed and " << skipped_alignments
			<< " alignments skipped (total: " << parsed_alignments + skipped_alignments
			<< ")" << endl;


		qint64 bam_file_parsing_runtime = timer.restart();
		out << "(BAM file parsing took " << get_time_string(bam_file_parsing_runtime)
			<< ")" << endl << endl;

		//Debug log:
		out << "summary: " << endl;
		out << "read counts\ttotal bases\tsoft-clipped bases\tinserted bases\tforward"
			<< "\treverse\tcorrect\tswitched" << endl;
		for(int chr = 0; chr < read_counts.size(); chr++) {
			long long rc = 0;
			long long tb = 0;
			long long scb = 0;
			long long ib = 0;
			long long fo = 0;
			long long ro = 0;
			long long co = 0;
			long long so = 0;
			for(int i= 0; i < read_counts[chr].size(); i++){
				rc += read_counts[chr][i];
				tb += total_bases[chr][i];
				scb += soft_clipped_bases[chr][i];
				ib += upper_insert_size_outliers[chr][i];
				fo += forward_orientated_reads[chr][i];
				ro += reverse_orientated_reads[chr][i];
				co += correct_orientated_reads[chr][i];
				so += switched_orientated_reads[chr][i];
			}
			out << rc << "\t" << tb << "\t" << scb << "\t" << ib << "\t" << fo << "\t" << ro << "\t" << co << "\t" << so << endl;
		}

		out << "writing output file ..." << endl;
		// write result to file
		QSharedPointer<QFile> outfile = Helper::openFileForWriting(getOutfile("out"), true);
		QTextStream outstream(outfile.data());

		// write header:
		outstream << "# " << getInfile("in") << endl;
		outstream << "#chr\tbinStart\tbinEnd\treadCounts\ttotalBases\tsoftclippedBases\tlowerInsertSizeOutliers\t"
				  << "upperInsertOutliers\tforwardOrientatedReads\treverseOrientatedReads\tcorrectOrientatedReads\tswitchedOrientatedReads"
				  << endl;
		// write data:
		for(int chr_id = 0; chr_id < chromosomes.size(); chr_id++){
			for(int bin_idx = 0; bin_idx < read_counts[chr_id].size(); bin_idx++){
				// write position
				outstream << reader.chromosome(chr_id).str() << "\t"
						  << (bin_idx * binSize) << "\t" << ((bin_idx + 1) * binSize) << "\t";
				// write data
				outstream << read_counts[chr_id][bin_idx] << "\t" << total_bases[chr_id][bin_idx] << "\t"
						  << soft_clipped_bases[chr_id][bin_idx] << "\t" << lower_insert_size_outliers[chr_id][bin_idx]
						  << "\t" << upper_insert_size_outliers[chr_id][bin_idx] << "\t"
						  << forward_orientated_reads[chr_id][bin_idx] << "\t" << reverse_orientated_reads[chr_id][bin_idx]
						  << "\t" << correct_orientated_reads[chr_id][bin_idx] << "\t"
						  << switched_orientated_reads[chr_id][bin_idx] << endl;

//				QString bed_line;
//				// write position
//				bed_line = reader.chromosome(chr_id).str() + "\t" + QString::number(bin_idx * binSize) + "\t" + QString::number((bin_idx + 1) * binSize) + "\t";
//				// write data
//				bed_line +=  QString::number(read_counts[chr_id][bin_idx]) + "\t" + QString::number(total_bases[chr_id][bin_idx]) + "\t"
//						  + QString::number(soft_clipped_bases[chr_id][bin_idx]) + "\t" + QString::number(lower_insert_size_outliers[chr_id][bin_idx])
//						  + "\t" + QString::number(upper_insert_size_outliers[chr_id][bin_idx]) + "\t"
//						  + QString::number(forward_orientated_reads[chr_id][bin_idx]) + "\t" + QString::number(reverse_orientated_reads[chr_id][bin_idx])
//						  + "\t" + QString::number(correct_orientated_reads[chr_id][bin_idx]) + "\t"
//						  + QString::number(switched_orientated_reads[chr_id][bin_idx]);

//				outstream << bed_line << endl;
			}
		}

		qint64 file_writing_runtime = timer.elapsed();
		out << "finished (runtime: " << get_time_string(file_writing_runtime) << ")" << endl << endl;

		out << "Total runtime: " << get_time_string(insert_size_cutoff_calculation_runtime
													+ bam_file_parsing_runtime
													+ file_writing_runtime)
			<< endl;

	}

private:
	int binSize;
	int minInsertSize;
	int maxInsertSize;
	int insertSizeDistributionSize;
	int minimumQualityScore;
};

#include "main.moc"

int main(int argc, char *argv[])
{
	ConcreteTool tool(argc, argv);
	return tool.execute();
}
