#ifndef REPORTWORKER_H
#define REPORTWORKER_H

#include <QObject>
#include <QString>
#include <QTextStream>
#include "VariantList.h"
#include "WorkerBase.h"
#include "ReportDialog.h"
#include "NGSD.h"
#include "FilterCascade.h"
#include "RtfDocument.h"

///Report generation worker.
class ReportWorker
		: public WorkerBase
{
	Q_OBJECT

public:
	///Constructor.
    ReportWorker(QString sample_name, QString file_bam, QString file_roi, const VariantList& variants, const FilterCascade& filters, ReportSettings settings, QStringList log_files, QString file_rep);
	virtual void process();

	///Returns the file to which the HTML report was written.
	QString getReportFile()
	{
		return file_rep_;
	}

	///writes a low-coverage report
	RtfTable writeCoverageReportTable(QString bam_file, QString roi_file, const BedFile& roi, const GeneSet& genes, int min_cov, NGSD& db, bool calculate_depth, QMap<QString, QString>* output=nullptr, bool gene_and_gap_details=true);
	RtfTable writeCoverageReportCCDS(QString bam_file, const GeneSet& genes, int min_cov, int extend, NGSD& db, QMap<QString, QString>* output=nullptr, bool gap_table=true, bool gene_details=true);

	///Returns if the pre-calcualed gaps for the given ROI.
	///If using the pre-calculated gaps file is not possible, @p message contains an error message.
	static BedFile precalculatedGaps(QString bam_file, const BedFile& roi, int min_cov, NGSD& db, QString& message);

	///Returns if given ROI file is the processing system target file corresponding to the BAM file.
	static bool isProcessingSystemTargetFile(QString bam_file, QString roi_file, NGSD& db);

	static void writeHtmlHeader(QTextStream& stream, QString sample_name);
	static void writeHtmlFooter(QTextStream& stream);
	static void validateAndCopyReport(QString from, QString to, bool put_to_archive, bool is_rtf);

	static QByteArray inheritance(const QByteArray& gene_info);

private:
	//input variables
	QString sample_name_;
	QString file_bam_;
	QString file_roi_;
	const VariantList& variants_;
    FilterCascade filters_;

	ReportSettings settings_;
	QStringList log_files_;
	//output variables
	QString file_rep_;

	//temporary variables
	BedFile roi_;
	GeneSet genes_;
	int var_count_;
	QMap<QString, QString> roi_stats_;

	//NGSD access
	NGSD db_;

	RtfSourceCode trans(const QByteArray& text) const;
	QList<QByteArray> formatCodingSplicing(const QList<VariantTranscript>& transcripts);
	QByteArray formatGenotype(const QByteArray& gender, const QByteArray& genotype, const Variant& variant);
	void writeHTML();
	void writeXML(QString outfile_name, QString report_document);


	void writeRtf();

	RtfDocument doc_;
};

#endif


