#include "ReportWorker.h"
#include "Log.h"
#include "Helper.h"
#include "Exceptions.h"
#include "Statistics.h"
#include "Settings.h"
#include "BedFile.h"
#include "ChromosomalIndex.h"
#include "VariantList.h"
#include "XmlHelper.h"
#include "NGSHelper.h"
#include "FilterCascade.h"
#include "GSvarHelper.h"

#include <QFile>
#include <QTextStream>
#include <QFileInfo>
#include <QCoreApplication>
#include <QXmlStreamWriter>
#include <QMessageBox>
#include <QDesktopServices>
#include <QApplication>



ReportWorker::ReportWorker(QString sample_name, QString file_bam, QString file_roi, const VariantList& variants, const FilterCascade& filters, ReportSettings settings, QStringList log_files, QString file_rep)
	: WorkerBase("Report generation")
	, sample_name_(sample_name)
	, file_bam_(file_bam)
	, file_roi_(file_roi)
	, variants_(variants)
	, filters_(filters)
	, settings_(settings)
	, log_files_(log_files)
	, file_rep_(file_rep)
	, var_count_(variants_.count())
{
}

void ReportWorker::process()
{
	//load ROI if given
	if (file_roi_!="")
	{
		roi_.load(file_roi_);
		roi_.merge();

		//determine variant count (inside target region)
		FilterResult filter_result(variants_.count());
		FilterRegions::apply(variants_, roi_, filter_result);
		var_count_ = filter_result.countPassing();

		//load gene list file
		genes_ = GeneSet::createFromFile(file_roi_.left(file_roi_.size()-4) + "_genes.txt");
	}

	roi_stats_.clear();

	writeRtf();
	//writeHTML();
}

QList<QByteArray> ReportWorker::formatCodingSplicing(const QList<VariantTranscript>& transcripts)
{
    const QMap<QByteArray, QByteArrayList>& preferred_transcripts = GSvarHelper::preferredTranscripts();

	QList<QByteArray> output;
	QList<QByteArray> output_pt;

	foreach(const VariantTranscript& trans, transcripts)
	{
		QByteArray line = trans.gene + ":" + trans.id + ":" + trans.hgvs_c + ":" + trans.hgvs_p;

		output.append(line);

        if (preferred_transcripts.value(trans.gene).contains(trans.id))
		{
			output_pt.append(line);
		}
	}

	//return only preferred transcripts if present
	if (output_pt.count()>0)
	{
		output = output_pt;
	}

	return output;
}

QByteArray ReportWorker::formatGenotype(const QByteArray& gender, const QByteArray& genotype, const Variant& variant)
{
	//correct only hom variants on gonosomes outside the PAR for males
	if (gender!="male") return genotype;
	if (genotype!="hom") return genotype;
	if (!variant.chr().isGonosome()) return genotype;
	if (NGSHelper::pseudoAutosomalRegion("hg19").overlapsWith(variant.chr(), variant.start(), variant.end())) return genotype;

	return "hemi";
}

RtfTable ReportWorker::writeCoverageReportTable(QString bam_file, QString roi_file, const BedFile& roi, const GeneSet& genes, int min_cov,  NGSD& db, bool calculate_depth, QMap<QString, QString>* output, bool gene_and_gap_details)
{
	//get target region coverages (from NGSD or calculate)
	QString avg_cov = "";
	QCCollection stats;
	if (isProcessingSystemTargetFile(bam_file, roi_file, db) || !calculate_depth)
	{
		try
		{
			QString processed_sample_id = db.processedSampleId(bam_file);
			stats = db.getQCData(processed_sample_id);
		}
		catch(...)
		{
		}
	}
	if (stats.count()==0)
	{
		Log::warn("Target region depth from NGSD cannot be used because ROI is not the processing system target region! Recalculating...");
		stats = Statistics::mapping(roi, bam_file);
	}
	for (int i=0; i<stats.count(); ++i)
	{
		if (stats[i].accession()=="QC:2000025") avg_cov = stats[i].toString();
	}

	RtfTable table;

	table.addRow(RtfTableRow({trans("Durchschnittliche Sequenziertiefe") + ":", avg_cov.toUtf8()},{3000,6636},RtfParagraph().setFontSize(16)));

	if (gene_and_gap_details)
	{
		//calculate low-coverage regions
		QString message;
		BedFile low_cov = precalculatedGaps(bam_file, roi, min_cov, db, message);
		if (!message.isEmpty())
		{
			Log::warn("Low-coverage statistics needs to be calculated. Pre-calulated gap file cannot be used because: " + message);
			low_cov = Statistics::lowCoverage(roi, bam_file, min_cov);
		}

		//annotate low-coverage regions with gene names
		for(int i=0; i<low_cov.count(); ++i)
		{
			BedLine& line = low_cov[i];
			GeneSet genes = db.genesOverlapping(line.chr(), line.start(), line.end(), 20); //extend by 20 to annotate splicing regions as well
			line.annotations().append(genes.join(", "));
		}

		//group by gene name
		QHash<QByteArray, BedFile> grouped;
		for (int i=0; i<low_cov.count(); ++i)
		{
			QList<QByteArray> genes = low_cov[i].annotations()[0].split(',');
			foreach(QByteArray gene, genes)
			{
				gene = gene.trimmed();

				//skip non-gene regions
				// - remains of VEGA database in old HaloPlex designs
				// - SNPs for sample identification
				if (gene=="") continue;

				grouped[gene].append(low_cov[i]);
			}
		}

		//output
		if (!genes.isEmpty())
		{
			QByteArrayList complete_genes;
			foreach(const QByteArray& gene, genes)
			{
				if (!grouped.contains(gene))
				{
					complete_genes << gene;
				}
			}
			table.addRow( RtfTableRow({trans("Komplett abgedeckte Gene") + ":", RtfText(complete_genes.join(", ")).setItalic(true).setFontSize(16).RtfCode()}, {3000,6636} , RtfParagraph().setFontSize(16) ) );
		}

		QByteArray gap_perc = QByteArray::number(100.0*low_cov.baseCount()/roi.baseCount(), 'f', 2);
		if (output!=nullptr) output->insert("gap_percentage", gap_perc);

		table.addRow(RtfTableRow({trans("Anteil Regionen mit Tiefe <") + " " + QByteArray::number(min_cov) + ":", gap_perc + "%"}, {3000,6636}, RtfParagraph().setFontSize(16)));
		if (!genes.isEmpty())
		{
			QList<RtfSourceCode> incomplete_genes;
			foreach(const QByteArray& gene, genes)
			{
				if (grouped.contains(gene))
				{
					incomplete_genes << gene + " " + RtfText(QByteArray::number(grouped[gene].baseCount())).setFontSize(14).RtfCode();
				}
			}
			table.addRow(RtfTableRow({trans("Fehlende Basen in nicht komplett abgedeckten Genen") + ":", incomplete_genes.join(", ")}, {3000,6636}, RtfParagraph().setFontSize(16)));
		}

		table.addRow(RtfTableRow(trans("Details Regionen mit Tiefe <") + QByteArray::number(min_cov),{9636}, RtfParagraph().setBold(true).setHorizontalAlignment("c")).setBackgroundColor(1));

		table.addRow(RtfTableRow({trans("Gen"), trans("Lücken"), trans("Chromosom"), trans("Koordinaten (hg19)")},{900,900,1200,6636}, RtfParagraph().setBold(true)).setHeader().setBackgroundColor(1));

		for (auto it=grouped.cbegin(); it!=grouped.cend(); ++it)
		{
			RtfTableRow row;

			const BedFile& gaps = it.value();
			QByteArray chr = gaps[0].chr().strNormalized(true);;
			QByteArrayList coords;
			for (int i=0; i<gaps.count(); ++i)
			{
				coords << QByteArray::number(gaps[i].start()) + "-" + QByteArray::number(gaps[i].end());
			}

			row.addCell(900, it.key(), RtfParagraph().setItalic(true).setFontSize(16));
			row.addCell(900, QByteArray::number(gaps.baseCount()), RtfParagraph().setFontSize(16));
			row.addCell(1200, chr, RtfParagraph().setFontSize(16));
			row.addCell(6636, coords.join(", "), RtfParagraph().setFontSize(16));

			table.addRow(row);
		}
	}


	table.setUniqueBorder(1,"brdrhair");

	return table;
}

RtfTable ReportWorker::writeCoverageReportCCDS(QString bam_file, const GeneSet& genes, int min_cov, int extend, NGSD& db, QMap<QString, QString>* output, bool gap_table, bool gene_details)
{
	QByteArray ext_string = (extend==0 ? "" : " ±" + QByteArray::number(extend) + " ");

	RtfTable table;

	if(gap_table) table.addRow(RtfTableRow(trans("Abdeckungsstatistik für CCDS") + " " + ext_string, 9636, RtfParagraph().setBold(true).setHorizontalAlignment("c")).setBackgroundColor(1));
	if(gap_table) table.addRow(RtfTableRow({trans("Gen"), trans("Transcript"), trans("Größe"), trans("Lücken"), trans("Chromosom"), trans("Koordinaten (hg19)") }, {900,1100,900,1000,1136,4600}, RtfParagraph().setBold(true)).setBackgroundColor(1));


	QMap<QByteArray, int> gap_count;
	long long bases_overall = 0;
	long long bases_sequenced = 0;
	GeneSet genes_noncoding;
	GeneSet genes_notranscript;
	foreach(const QByteArray& gene, genes)
	{
		int gene_id = db.geneToApprovedID(gene);

		//approved gene symbol
		QByteArray symbol = db.geneSymbol(gene_id);

		//longest coding transcript
		Transcript transcript = db.longestCodingTranscript(gene_id, Transcript::CCDS, true);
		if (!transcript.isValid())
		{
			transcript = db.longestCodingTranscript(gene_id, Transcript::CCDS, true, true);
			if (!transcript.isValid() || transcript.regions().baseCount()==0)
			{
				genes_notranscript.insert(gene);

				if(gap_table) table.addRow(RtfTableRow({symbol, "n/a", "n/a", "n/a", "n/a", "n/a"},  {900,1100,900,1000,1136,4600}, RtfParagraph().setFontSize(16))) ;
				continue;
			}
			else
			{
				genes_noncoding.insert(gene);
			}
		}

		//gaps
		QString message;
		BedFile roi = transcript.regions();
		if (extend>0)
		{
			roi.extend(extend);
			roi.merge();
		}
		BedFile gaps = precalculatedGaps(bam_file, roi, min_cov, db, message);
		if (!message.isEmpty())
		{
			Log::warn("Low-coverage statistics for transcript " + transcript.name() + " needs to be calculated. Pre-calulated gap file cannot be used because: " + message);
			gaps = Statistics::lowCoverage(roi, bam_file, min_cov);
		}

		long long bases_transcipt = roi.baseCount();
		long long bases_gaps = gaps.baseCount();
		QByteArrayList coords;
		for (int i=0; i<gaps.count(); ++i)
		{
			coords << QByteArray::number(gaps[i].start()) + "-" + QByteArray::number(gaps[i].end());
		}

		if(gap_table) table.addRow(RtfTableRow({symbol, transcript.name(), QByteArray::number(bases_transcipt), QByteArray::number(bases_gaps), roi[0].chr().strNormalized(true), coords.join(", ")},{900,1100,900,1000,1136,4600},RtfParagraph().setFontSize(16)));
		gap_count[symbol] += bases_gaps;
		bases_overall += bases_transcipt;
		bases_sequenced += bases_transcipt - bases_gaps;
	}

	//show warning if non-coding transcripts had to be used
	if (!genes_noncoding.isEmpty())
	{
		table.addRow(RtfTableRow("Warning: Using the longest *non-coding* transcript for genes " + genes_noncoding.join(", ") + " (no coding transcripts for GRCh37 defined)", 9636, RtfParagraph().setFontSize(16)));
	}
	if (!genes_notranscript.isEmpty())
	{
		table.addRow(RtfTableRow("Warning: No transcript defined for genes " + genes_notranscript.join(", "), 9636, RtfParagraph().setFontSize(16)));
	}

	table.setUniqueBorder(1,"brdrhair");

	//overall statistics
	table.addRow(RtfTableRow(" ",9636));

	table.addRow(RtfTableRow({"CCDS" + ext_string + trans("gesamt") + ":", QByteArray::number(bases_overall)}, {3000, 6636}, RtfParagraph().setFontSize(16) ) );
	table.addRow(RtfTableRow({"CCDS" + ext_string + trans("mit Tiefe") + " < " + QByteArray::number(min_cov) + ":", QByteArray::number(bases_sequenced) + " (" + QByteArray::number(100.0 * bases_sequenced / bases_overall, 'f', 2) + "%)"}, {3000, 6636}, RtfParagraph().setFontSize(16) ) );

	long long gaps = bases_overall - bases_sequenced;

	table.addRow(RtfTableRow({"CCDS" + ext_string + trans("mit Tiefe") + " < " + QByteArray::number(min_cov) + ":", QByteArray::number(gaps) + " (" + QByteArray::number(100.0 * gaps / bases_overall, 'f', 2) + "%)"}, {3000, 6636}, RtfParagraph().setFontSize(16) ) );

	//gene statistics
	if (gene_details)
	{
		QList<RtfSourceCode> genes_complete;
		QByteArrayList genes_incomplete;
		for (auto it = gap_count.cbegin(); it!=gap_count.cend(); ++it)
		{
			if (it.value()==0)
			{
				genes_complete << it.key();
			}
			else
			{
				genes_incomplete << it.key() + " " + RtfText(QByteArray::number(it.value())).setFontSize(14).RtfCode();
			}
		}
		table.addRow(RtfTableRow({trans("Komplett abgedeckte Gene") + ":", genes_complete.join(", ")}, {3000, 6636}, RtfParagraph().setFontSize(16)));

		table.addRow(RtfTableRow({trans("Fehlende Basen in nicht komplett abgedeckten Genen") + ":", genes_incomplete.join(", ")}, {3000, 6636}, RtfParagraph().setFontSize(16)));
	}

	if (output!=nullptr) output->insert("ccds_sequenced", QString::number(bases_sequenced));


	return table;
}

BedFile ReportWorker::precalculatedGaps(QString bam_file, const BedFile& roi, int min_cov, NGSD& db, QString& message)
{
	message.clear();

	//check depth cutoff
	if (min_cov!=20)
	{
		message = "Depth cutoff is not 20!";
		return BedFile();
	}

	//find low-coverage file
	QString dir = QFileInfo(bam_file).absolutePath();
	QStringList low_cov_files = Helper::findFiles(dir, "*_lowcov.bed", false);
	if(low_cov_files.count()!=1)
	{
		message = "Low-coverage file does not exist in " + dir;
		return BedFile();
	}
	QString low_cov_file = low_cov_files[0];

	//load low-coverage file
	BedFile gaps;
	gaps.load(low_cov_file);

	//For WGS there is nothing more to check
	QString processed_sample_id = db.processedSampleId(bam_file);
	ProcessingSystemData system_data = db.getProcessingSystemData(processed_sample_id, true);
	if (system_data.type=="WGS")
	{
		return gaps;
	}

	//extract processing system ROI statistics
	int regions = -1;
	long long bases = -1;
	foreach(QString line, gaps.headers())
	{
		if (line.startsWith("#ROI bases: "))
		{
			bool ok = true;
			bases = line.mid(12).trimmed().toLongLong(&ok);
			if (!ok) bases = -1;
		}
		if (line.startsWith("#ROI regions: "))
		{
			bool ok = true;
			regions = line.mid(14).trimmed().toInt(&ok);
			if (!ok) regions = -1;
		}
	}
	if (regions<0 || bases<0)
	{
		message = "Low-coverage file header does not contain target region statistics: " + low_cov_file;
		return BedFile();
	}

	//compare statistics to current processing system
	if (system_data.target_file=="")
	{
		message = "Processing system target file not defined in NGSD!";
		return BedFile();
	}
	BedFile sys;
	sys.load(system_data.target_file);
	sys.merge();
	if (sys.count()!=regions || sys.baseCount()!=bases)
	{
		message = "Low-coverage file is outdated. It does not match processing system target region: " + low_cov_file;
		return BedFile();
	}

	//calculate gaps inside target region
	gaps.intersect(roi);

	//add target region bases not covered by processing system target file
	BedFile uncovered(roi);
	uncovered.subtract(sys);
	gaps.add(uncovered);
	gaps.merge();

	return gaps;
}

bool ReportWorker::isProcessingSystemTargetFile(QString bam_file, QString roi_file, NGSD& db)
{
	ProcessingSystemData system_data = db.getProcessingSystemData(db.processedSampleId(bam_file), true);

	return Helper::canonicalPath(system_data.target_file) == Helper::canonicalPath(roi_file);
}

void ReportWorker::writeHtmlHeader(QTextStream& stream, QString sample_name)
{
	stream << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">" << endl;
	stream << "<html xmlns=\"http://www.w3.org/1999/xhtml\">" << endl;
	stream << "	<head>" << endl;
	stream << "	   <title>Report " << sample_name << "</title>" << endl;
	stream << "	   <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />" << endl;
	stream << "	   <style type=\"text/css\">" << endl;
	stream << "		<!--" << endl;
	stream << "body" << endl;
	stream << "{" << endl;
	stream << "	font-family: sans-serif;" << endl;
	stream << "	font-size: 70%;" << endl;
	stream << "}" << endl;
	stream << "table" << endl;
	stream << "{" << endl;
	stream << "	border-collapse: collapse;" << endl;
	stream << "	border: 1px solid black;" << endl;
	stream << "	width: 100%;" << endl;
	stream << "}" << endl;
	stream << "th, td" << endl;
	stream << "{" << endl;
	stream << "	border: 1px solid black;" << endl;
	stream << "	font-size: 100%;" << endl;
	stream << "	text-align: left;" << endl;
	stream << "}" << endl;
	stream << "p" << endl;
	stream << "{" << endl;
	stream << " margin-bottom: 0cm;" << endl;
	stream << "}" << endl;
	stream << "		-->" << endl;
	stream << "	   </style>" << endl;
	stream << "	</head>" << endl;
	stream << "	<body>" << endl;
}

void ReportWorker::writeHtmlFooter(QTextStream& stream)
{
	stream << "	</body>" << endl;
	stream << "</html>" << endl;
}


void ReportWorker::writeRtf()
{
	doc_.addColor(217,217,217);
	doc_.addColor(255,0,0);

	//get data from database
	QString sample_id = db_.sampleId(sample_name_);
	SampleData sample_data = db_.getSampleData(sample_id);
	QString processed_sample_id = db_.processedSampleId(sample_name_);
	ProcessedSampleData processed_sample_data = db_.getProcessedSampleData(processed_sample_id);
	ProcessingSystemData system_data = db_.getProcessingSystemData(processed_sample_id, true);


	doc_.addPart(RtfParagraph("Technischer Report zur bioinformatischen Analyse").setFontSize(18).setBold(true).RtfCode());

	doc_.addPart(RtfParagraph(RtfText(trans("Probe") + ": " + sample_name_.toUtf8()).setBold(true).RtfCode() + " (" + sample_data.name_external.toUtf8() + ")").RtfCode());

	QList<RtfSourceCode> text;

	text << trans("Geschlecht") + ": " + processed_sample_data.gender.toUtf8();
	text << trans("Prozessierungssystem") + ": " + processed_sample_data.processing_system.toUtf8();
	text << trans("Referenzgenom") + ": " + system_data.genome.toUtf8();
	text << trans("Datum") + ": " + QDate::currentDate().toString("dd.MM.yyyy").toUtf8();
	text << trans("Benutzer") + ": " + Helper::userName().toUtf8();
	text << trans("Analysepipeline") + ": " + variants_.getPipeline().toUtf8();
	text << trans("Auswertungssoftware") + ": " + QCoreApplication::applicationName().toUtf8() + " " + QCoreApplication::applicationVersion().toUtf8();
	text << trans("KASP-Ergebnis") + ": " + db_.getQCData(processed_sample_id).value("kasp").asString().toUtf8();

	doc_.addPart(RtfParagraph(RtfText(text).setBold(true).RtfCode()).RtfCode());
	text.clear();


	//Phenotype information
	doc_.addPart(RtfParagraph(trans("Phänotyp")).setFontSize(18).setBold(true).RtfCode());

	QList<SampleDiseaseInfo> info = db_.getSampleDiseaseInfo(sample_id, "ICD10 code");
	foreach(const SampleDiseaseInfo& entry, info)
	{
		text << "ICD10: " + entry.disease_info.toUtf8();
	}

	info = db_.getSampleDiseaseInfo(sample_id, "HPO term id");
	foreach(const SampleDiseaseInfo& entry, info)
	{
		text << "HPO: " << entry.disease_info.toUtf8() + " (" + db_.phenotypeByAccession(entry.disease_info.toLatin1(), false).name() + ")";
	}
	info = db_.getSampleDiseaseInfo(sample_id, "OMIM disease/phenotype identifier");
	foreach(const SampleDiseaseInfo& entry, info)
	{
		text << "OMIM: " + entry.disease_info.toUtf8();
	}
	info = db_.getSampleDiseaseInfo(sample_id, "Orpha number");
	foreach(const SampleDiseaseInfo& entry, info)
	{
		text << "Orphanet: " << entry.disease_info.toUtf8();
	}
	doc_.addPart(RtfParagraph(RtfText(text).RtfCode()).RtfCode());
	text.clear();



	//Target region statistics
	if (file_roi_!="")
	{
		doc_.addBlankLine();
		doc_.addPart(RtfParagraph(trans("Zielregion")).setFontSize(18).setBold(true).RtfCode());


		doc_.addPart(RtfParagraph(trans("Die Zielregion umfasst mindestens die CCDS (\"consensus coding sequence\") unten genannter Gene ±20 Basen flankierender intronischer Sequenz, kann aber auch zusätzliche Exons und/oder flankierende Basen beinhalten.")).setFontSize(14).RtfCode());
		doc_.addPart(RtfParagraph(RtfText(trans("Name:")).setBold(true).setFontSize(14).RtfCode() + " " + QFileInfo(file_roi_).fileName().replace(".bed", "").toUtf8()).setFontSize(14).RtfCode());
		if (!genes_.isEmpty())
		{
			 doc_.addPart(RtfParagraph(RtfText(trans("Ausgewertete Gene") + " (" + QByteArray::number(genes_.count()) + "):").setBold(true).setFontSize(14).RtfCode() + " " + genes_.join(", ")).setFontSize(14).RtfCode());
		}
	}


	//get column indices
	int i_genotype = variants_.getSampleHeader().infoByStatus(true).column_index;
	int i_gene = variants_.annotationIndexByName("gene", true, true);
	int i_co_sp = variants_.annotationIndexByName("coding_and_splicing", true, true);
	int i_omim = variants_.annotationIndexByName("OMIM", true, true);
	int i_class = variants_.annotationIndexByName("classification", true, true);
	int i_comment = variants_.annotationIndexByName("comment", true, true);
	int i_kg = variants_.annotationIndexByName("1000G", true, true);
	int i_gnomad = variants_.annotationIndexByName("gnomAD", true, true);

	//output: applied filters
	doc_.addBlankLine();
	doc_.addPart(RtfParagraph(trans("Filterkriterien")).setBold(true).RtfCode());

	doc_.addPart(RtfParagraph(trans("Gefundene Varianten in Zielregion gesamt") + ": " + QByteArray::number(var_count_)).RtfCode());
	doc_.addPart(RtfParagraph(trans("Anzahl Varianten ausgewählt für Report") + ": " + QByteArray::number(settings_.report_config.variantIndices(VariantType::SNVS_INDELS, true, settings_.report_type).count())).RtfCode());

	for(int i=0; i<filters_.count(); ++i)
	{
		text << "- " + filters_[i]->toText().replace(",",", ").toUtf8();
	}
	doc_.addPart(RtfParagraph(RtfText(text).RtfCode()).setIndent(50,0,0).RtfCode());

	doc_.addPart(RtfParagraph("").RtfCode());
	doc_.addPart(RtfParagraph(trans("Varianten nach klinischer Interpretation im Kontext der Fragestellung")).setBold(true).setFontSize(18).RtfCode());

	RtfTable var_table;

	var_table.addRow(RtfTableRow({trans("Gen"), trans("Variante"),  trans("Genotyp"), trans("Details"), trans("Typ"), trans("Klasse"), trans("Vererbung"), "1000g", "gnomAD"},{900,2100,800,2400,500,650,887,600,800},RtfParagraph().setFontSize(16).setBold(true).setHorizontalAlignment("c")).setBackgroundColor(1).setHeader());


	for(const auto& var_conf : settings_.report_config.variantConfig())
	{
		if (var_conf.variant_type!=VariantType::SNVS_INDELS) continue;
		if (!var_conf.showInReport()) continue;
		if (var_conf.report_type!=settings_.report_type) continue;

		const Variant& variant = variants_[var_conf.variant_index];
		QByteArray genes = variant.annotations()[i_gene];

		RtfTableRow row;
		row.addCell(900, genes, RtfParagraph().setItalic(true).setFontSize(16));
		row.addCell(2100, variant.chr().str() + ":" + QByteArray::number(variant.start()) + " " + variant.ref()  + " > " + variant.obs(), RtfParagraph().setFontSize(16));
		row.addCell(800, formatGenotype(processed_sample_data.gender.toLatin1(), variant.annotations().at(i_genotype), variant), RtfParagraph().setFontSize(16));
		row.addCell(formatCodingSplicing(variant.transcriptAnnotations(i_co_sp)), 2400, RtfParagraph().setFontSize(16));


		QByteArrayList type_strings;
		type_strings << var_conf.report_type.toUtf8().replace("report: ", "");
		if (var_conf.de_novo) type_strings << "de-novo";
		if (var_conf.mosaic) type_strings << "mosaic";
		if (var_conf.comp_het) type_strings << "comp-het";

		row.addCell(500, type_strings.join(","), RtfParagraph().setFontSize(16));
		row.addCell(650, variant.annotations().at(i_class), RtfParagraph().setFontSize(16));
		row.addCell(887, var_conf.inheritance.toUtf8(), RtfParagraph().setFontSize(16));

		QByteArray freq = variant.annotations().at(i_kg).trimmed();
		row.addCell(600, freq.isEmpty() ? "n/a" : freq, RtfParagraph().setFontSize(16));
		freq = variant.annotations().at(i_gnomad).trimmed();
		row.addCell(800, freq.isEmpty() ? "n/a" : freq, RtfParagraph().setFontSize(16));

		var_table.addRow(row);


		//OMIM and comment line
		QString omim = variant.annotations()[i_omim];
		QByteArray comment = variant.annotations()[i_comment];
		if (comment!="" || omim!="")
		{
			QList<RtfSourceCode> parts;
			if (comment!="") parts << RtfText("NGSD: " + comment).highlight(2).setFontSize(16).RtfCode();
			if (omim!="")
			{
				QStringList omim_parts = omim.append(" ").split("]; ");
				for(QString omim_part : omim_parts)
				{
					if (omim_part.count()<10) continue;
					omim = "OMIM ID: " + omim_part.left(6).toUtf8() + " Details: " + omim_part.mid(8).toUtf8();
				}

				parts << omim.replace(",",", ").toUtf8();
			}

			var_table.addRow(RtfTableRow(parts.join("\\line\\n"),{9637},RtfParagraph().setFontSize(16)));
		}

	}

	var_table.setUniqueBorder(1,"brdrhair");

	doc_.addPart(var_table.RtfCode());

	doc_.addPart(RtfParagraph("").RtfCode());

	doc_.addPart(RtfParagraph(trans("Für Informationen zur Klassifizierung von Varianten, siehe allgemeine Zusatzinformationen.")).RtfCode());
	doc_.addPart(RtfParagraph(trans("Teilweise können bei Varianten unklarer Signifikanz (Klasse 3) - in Abhängigkeit von der Art der genetischen Veränderung, der Familienanamnese und der Klinik des/der Patienten - weiterführende Untersuchungen eine änderung der Klassifizierung bewirken. Bei konkreten differentialdiagnostischen Hinweisen auf eine entsprechende Erkrankung ist eine humangenetische Mitbeurteilung erforderlich, zur Beurteilung ob erweiterte genetische Untersuchungen zielführend wären.")).setHorizontalAlignment("j").RtfCode());

	//classification explanation
	if(settings_.show_class_details)
	{
		text.clear();
		text << RtfText(trans("Klassifikation von Varianten") + ":").setBold(true).RtfCode();
		text << trans("Die Klassifikation der Varianten erfolgt in Anlehnung an die Publikation von Plon et al. (Hum Mutat 2008)");
		text << RtfText(trans("Klasse 5: Eindeutig pathogene Veränderung / Mutation") + ":").setBold(true).RtfCode() + " " + trans("Veränderung, die bereits in der Fachliteratur mit ausreichender Evidenz als krankheitsverursachend bezogen auf das vorliegende Krankheitsbild beschrieben wurde sowie als pathogen zu wertende Mutationstypen (i.d.R. Frameshift- bzw. Stoppmutationen).");
		text << RtfText(trans("Klasse 4: Wahrscheinlich pathogene Veränderung") + ":").setBold(true).RtfCode() + " " + trans("DNA-Veränderung, die aufgrund ihrer Eigenschaften als sehr wahrscheinlich krankheitsverursachend zu werten ist.");
		text << RtfText(trans("Klasse 3: Variante unklarer Signifikanz (VUS) - Unklare Pathogenität") + ":").setBold(true).RtfCode() + " " + trans("Variante, bei der es unklar ist, ob eine krankheitsverursachende Wirkung besteht. Diese Varianten werden tabellarisch im technischen Report mitgeteilt.");
		text << RtfText(trans("Klasse 2: Sehr wahrscheinlich benigne Veränderungen") + ":").setBold(true).RtfCode() + " " + trans("Aufgrund der Häufigkeit in der Allgemeinbevölkerung oder der Lokalisation bzw. aufgrund von Angaben in der Literatur sehr wahrscheinlich benigne. Werden nicht mitgeteilt, können aber erfragt werden.");
		text << RtfText(trans("Klasse 1: Benigne Veränderungen") + ":").setBold(true).RtfCode() + " " + trans("Werden nicht mitgeteilt, können aber erfragt werden.");

		doc_.addPart(RtfParagraph(text.join("\\line\n")).setFontSize(16).RtfCode());
	}

	//low-coverage analysis
	if(settings_.show_coverage_details && file_bam_ != "")
	{
		doc_.addPart(RtfParagraph("").RtfCode());
		doc_.addPart(RtfParagraph(trans("Abdeckungsstatistik")).setBold(true).RtfCode());
		doc_.addPart(writeCoverageReportTable(file_bam_, file_roi_, roi_, genes_, settings_.min_depth, db_, settings_.recalculate_avg_depth, &roi_stats_, settings_.roi_low_cov).RtfCode());


		doc_.addBlankLine();

		doc_.addPart(RtfParagraph().RtfCode());

		writeCoverageReportCCDS(file_bam_, genes_, settings_.min_depth, 0, db_, &roi_stats_, false, false);
		doc_.addPart(writeCoverageReportCCDS(file_bam_, genes_, settings_.min_depth, 5, db_, nullptr, true, true).RtfCode());
	}


	//OMIM table
	if (settings_.show_omim_table)
	{
		//prepare queries
		SqlQuery q_genes = db_.getQuery();
		q_genes.prepare("SELECT id, mim FROM omim_gene WHERE gene=:1");

		RtfTable table;
		table.addRow(RtfTableRow(trans("OMIM Gene und Phänotypen"), 9636, RtfParagraph().setBold(true).setHorizontalAlignment("c")).setBackgroundColor(1));
		table.addRow(RtfTableRow({trans("Gen"), trans("OMIM Gen MIM"), trans("OMIM Phänotypen")}, {1200, 1200, 7236}, RtfParagraph().setBold(true)).setBackgroundColor(1));

		foreach(QByteArray gene, genes_)
		{
			//approved gene symbol
			QByteArray gene_approved = db_.geneToApproved(gene, true);

			//generate table rows
			q_genes.bindValue(0, gene_approved);
			q_genes.exec();
			while (q_genes.next())
			{
				QByteArray gene_id = q_genes.value(0).toByteArray();
				QByteArray gene_mim = q_genes.value(1).toByteArray();
				RtfSourceCode phenotype = db_.getValues("SELECT phenotype FROM omim_phenotype WHERE omim_gene_id=" + gene_id).join("\\line\n").toUtf8();
				table.addRow(RtfTableRow({gene, gene_mim,phenotype},{1200,1200,7236},RtfParagraph().setFontSize(16)));
			}
		}
		table.setUniqueBorder(1,"brdrhair");

		doc_.addBlankLine();
		doc_.addPart(table.RtfCode());
	}


	doc_.save("C:\\Users\\ahgscha1\\Desktop\\rtf_germline_test.rtf");
}

void ReportWorker::writeHTML()
{
	QString temp_filename = Helper::tempFileName(".html");
	QSharedPointer<QFile> outfile = Helper::openFileForWriting(temp_filename);
	QTextStream stream(outfile.data());
	writeHtmlHeader(stream, sample_name_);

	//get trio data
	bool is_trio = variants_.type() == GERMLINE_TRIO;
	SampleInfo info_father;
	SampleInfo info_mother;
	if (is_trio)
	{
		info_father = variants_.getSampleHeader().infoByStatus(false, "male");
		info_mother = variants_.getSampleHeader().infoByStatus(false, "female");
	}

	//get data from database
	QString sample_id = db_.sampleId(sample_name_);
	SampleData sample_data = db_.getSampleData(sample_id);
	QString processed_sample_id = db_.processedSampleId(sample_name_);
	ProcessedSampleData processed_sample_data = db_.getProcessedSampleData(processed_sample_id);
	ProcessingSystemData system_data = db_.getProcessingSystemData(processed_sample_id, true);

	//report header (meta information)
	stream << "<h4>" << trans("Technischer Report zur bioinformatischen Analyse") << "</h4>" << endl;

	stream << "<p>" << endl;
	stream << "<b>" << trans("Probe") << ": " << sample_name_ << "</b> (" << sample_data.name_external << ")" << endl;
	if (is_trio)
	{
		stream << "<br />" << endl;
		stream << "<br />" << trans("Vater") << ": "  << info_father.id << endl;
		stream << "<br />" << trans("Mutter") << ": "  << info_mother.id << endl;
	}
	stream << "<br />" << endl;
	stream << "<br />" << trans("Geschlecht") << ": " << processed_sample_data.gender << endl;
	stream << "<br />" << trans("Prozessierungssystem") << ": " << processed_sample_data.processing_system << endl;
	stream << "<br />" << trans("Referenzgenom") << ": " << system_data.genome << endl;
	stream << "<br />" << trans("Datum") << ": " << QDate::currentDate().toString("dd.MM.yyyy") << endl;
	stream << "<br />" << trans("Benutzer") << ": " << Helper::userName() << endl;
	stream << "<br />" << trans("Analysepipeline") << ": "  << variants_.getPipeline() << endl;
	stream << "<br />" << trans("Auswertungssoftware") << ": "  << QCoreApplication::applicationName() << " " << QCoreApplication::applicationVersion() << endl;
	stream << "<br />" << trans("KASP-Ergebnis") << ": " << db_.getQCData(processed_sample_id).value("kasp").asString() << endl;
	stream << "</p>" << endl;

	///Phenotype information
	stream << "<p><b>" << trans("Phänotyp") << "</b>" << endl;
	QList<SampleDiseaseInfo> info = db_.getSampleDiseaseInfo(sample_id, "ICD10 code");
	foreach(const SampleDiseaseInfo& entry, info)
	{
		stream << "<br />ICD10: " << entry.disease_info << endl;
	}
	info = db_.getSampleDiseaseInfo(sample_id, "HPO term id");
	foreach(const SampleDiseaseInfo& entry, info)
	{
		stream << "<br />HPO: " << entry.disease_info << " (" << db_.phenotypeByAccession(entry.disease_info.toLatin1(), false).name() << ")" << endl;
	}
	info = db_.getSampleDiseaseInfo(sample_id, "OMIM disease/phenotype identifier");
	foreach(const SampleDiseaseInfo& entry, info)
	{
		stream << "<br />OMIM: " << entry.disease_info << endl;
	}
	info = db_.getSampleDiseaseInfo(sample_id, "Orpha number");
	foreach(const SampleDiseaseInfo& entry, info)
	{
		stream << "<br />Orphanet: " << entry.disease_info << endl;
	}
	stream << "</p>" << endl;

	///Target region statistics
	if (file_roi_!="")
	{
		stream << "<p><b>" << trans("Zielregion") << "</b>" << endl;
		stream << "<br /><span style=\"font-size: 80%;\">" << trans("Die Zielregion umfasst mindestens die CCDS (\"consensus coding sequence\") unten genannter Gene &plusmn;20 Basen flankierender intronischer Sequenz, kann aber auch zusätzliche Exons und/oder flankierende Basen beinhalten.") << endl;
		stream << "<br />" << trans("Name") << ": " << QFileInfo(file_roi_).fileName().replace(".bed", "") << endl;
		if (!genes_.isEmpty())
		{
			stream << "<br />" << trans("Ausgewertete Gene") << " (" << QString::number(genes_.count()) << "): " << genes_.join(", ") << endl;
		}
		stream << "</span></p>" << endl;
	}

	//get column indices
	int i_genotype = variants_.getSampleHeader().infoByStatus(true).column_index;
	int i_gene = variants_.annotationIndexByName("gene", true, true);
	int i_co_sp = variants_.annotationIndexByName("coding_and_splicing", true, true);
	int i_omim = variants_.annotationIndexByName("OMIM", true, true);
	int i_class = variants_.annotationIndexByName("classification", true, true);
	int i_comment = variants_.annotationIndexByName("comment", true, true);
	int i_kg = variants_.annotationIndexByName("1000G", true, true);
	int i_gnomad = variants_.annotationIndexByName("gnomAD", true, true);

	//output: applied filters
	stream << "<p><b>" << trans("Filterkriterien") << " " << "</b>" << endl;
	stream << "<br />" << trans("Gefundene Varianten in Zielregion gesamt") << ": " << var_count_ << endl;
	stream << "<br />" << trans("Anzahl Varianten ausgewählt für Report") << ": " << settings_.report_config.variantIndices(VariantType::SNVS_INDELS, true, settings_.report_type).count() << endl;
	for(int i=0; i<filters_.count(); ++i)
	{
		stream << "<br />&nbsp;&nbsp;&nbsp;&nbsp;- " << filters_[i]->toText() << endl;
	}
	stream << "</p>" << endl;

	//output: selected variants
	stream << "<p><b>" << trans("Varianten nach klinischer Interpretation im Kontext der Fragestellung") << "</b>" << endl;
	stream << "</p>" << endl;
	stream << "<table>" << endl;
	stream << "<tr><td><b>" << trans("Gen") << "</b></td><td><b>" << trans("Variante") << "</b></td><td><b>" << trans("Genotyp") << "</b></td>";
	if (is_trio)
	{
		stream << "<td><b>" << trans("Vater") << "</b></td>";
		stream << "<td><b>" << trans("Mutter") << "</b></td>";
	}
	stream << "<td><b>" << trans("Details") << "</b></td><td><b>" << trans("Klasse") << "</b></td><td><b>" << trans("Vererbung") << "</b></td><td><b>1000g</b></td><td><b>gnomAD</b></td></tr>" << endl;

	foreach(const ReportVariantConfiguration& var_conf, settings_.report_config.variantConfig())
	{
		if (var_conf.variant_type!=VariantType::SNVS_INDELS) continue;
		if (!var_conf.showInReport()) continue;
		if (var_conf.report_type!=settings_.report_type) continue;

		const Variant& variant = variants_[var_conf.variant_index];
		QByteArray genes = variant.annotations()[i_gene];
		stream << "<tr>" << endl;
		stream << "<td>" << genes << "</td>" << endl;
		stream << "<td>" << endl;
		stream  << variant.chr().str() << ":" << variant.start() << "&nbsp;" << variant.ref() << "&nbsp;>&nbsp;" << variant.obs() << "</td>";
		QStringList geno_info;
		geno_info << formatGenotype(processed_sample_data.gender.toLatin1(), variant.annotations().at(i_genotype), variant);
		if (var_conf.de_novo) geno_info << "de-novo";
		if (var_conf.mosaic) geno_info << "mosaic";
		if (var_conf.comp_het) geno_info << "comp-het";
		stream << "<td>" << geno_info.join(", ") << "</td>" << endl;
		if (is_trio)
		{
			stream << "<td>" << formatGenotype("male", variant.annotations().at(info_father.column_index), variant) << "</td>";
			stream << "<td>" << formatGenotype("female", variant.annotations().at(info_mother.column_index), variant) << "</td>";
		}
		//stream << "<td>" << formatCodingSplicing(variant.transcriptAnnotations(i_co_sp)) << "</td>" << endl;
		stream << "<td>" << variant.annotations().at(i_class) << "</td>" << endl;
		stream << "<td>" << var_conf.inheritance << "</td>" << endl;
		QByteArray freq = variant.annotations().at(i_kg).trimmed();
		stream << "<td>" << (freq.isEmpty() ? "n/a" : freq) << "</td>" << endl;
		freq = variant.annotations().at(i_gnomad).trimmed();
		stream << "<td>" << (freq.isEmpty() ? "n/a" : freq) << "</td>" << endl;
		stream << "</tr>" << endl;

		//OMIM and comment line
		QString omim = variant.annotations()[i_omim];
		QString comment = variant.annotations()[i_comment];
		if (comment!="" || omim!="")
		{
			QStringList parts;
			if (comment!="") parts << "<span style=\"background-color: #FF0000\">NGSD: " + comment + "</span>";
			if (omim!="")
			{
				QStringList omim_parts = omim.append(" ").split("]; ");
				foreach(QString omim_part, omim_parts)
				{
					if (omim_part.count()<10) continue;
					omim = "OMIM ID: " + omim_part.left(6) + " Details: " + omim_part.mid(8);
				}

				parts << omim;
			}
			stream << "<tr><td colspan=\"9\">" << parts.join("<br />") << "</td></tr>" << endl;
		}
	}
	stream << "</table>" << endl;

	stream << "<p>" << trans("Für Informationen zur Klassifizierung von Varianten, siehe allgemeine Zusatzinformationen.") << endl;
	stream << "</p>" << endl;

	stream << "<p>" << trans("Teilweise können bei Varianten unklarer Signifikanz (Klasse 3) -  in Abhängigkeit von der Art der genetischen Veränderung, der Familienanamnese und der Klinik des/der Patienten - weiterführende Untersuchungen eine änderung der Klassifizierung bewirken. Bei konkreten differentialdiagnostischen Hinweisen auf eine entsprechende Erkrankung ist eine humangenetische Mitbeurteilung erforderlich, zur Beurteilung ob erweiterte genetische Untersuchungen zielführend wären.") << endl;
	stream << "</p>" << endl;

	///classification explaination
	if (settings_.show_class_details)
	{
		stream << "<p><b>" << trans("Klassifikation von Varianten") << ":</b>" << endl;
		stream << "<br />" << trans("Die Klassifikation der Varianten erfolgt in Anlehnung an die Publikation von Plon et al. (Hum Mutat 2008)") << endl;
		stream << "<br /><b>" << trans("Klasse 5: Eindeutig pathogene Veränderung / Mutation") << ":</b> " << trans("Veränderung, die bereits in der Fachliteratur mit ausreichender Evidenz als krankheitsverursachend bezogen auf das vorliegende Krankheitsbild beschrieben wurde sowie als pathogen zu wertende Mutationstypen (i.d.R. Frameshift- bzw. Stoppmutationen).") << endl;
		stream << "<br /><b>" << trans("Klasse 4: Wahrscheinlich pathogene Veränderung") << ":</b> " << trans("DNA-Veränderung, die aufgrund ihrer Eigenschaften als sehr wahrscheinlich krankheitsverursachend zu werten ist.") << endl;
		stream << "<br /><b>" << trans("Klasse 3: Variante unklarer Signifikanz (VUS) - Unklare Pathogenität") << ":</b> " << trans("Variante, bei der es unklar ist, ob eine krankheitsverursachende Wirkung besteht. Diese Varianten werden tabellarisch im technischen Report mitgeteilt.") << endl;
		stream << "<br /><b>" << trans("Klasse 2: Sehr wahrscheinlich benigne Veränderungen") << ":</b> " << trans("Aufgrund der Häufigkeit in der Allgemeinbevölkerung oder der Lokalisation bzw. aufgrund von Angaben in der Literatur sehr wahrscheinlich benigne. Werden nicht mitgeteilt, können aber erfragt werden.") << endl;
		stream << "<br /><b>" << trans("Klasse 1: Benigne Veränderungen") << ":</b> " << trans("Werden nicht mitgeteilt, können aber erfragt werden.") << endl;
		stream << "</p>" << endl;
	}

	///low-coverage analysis
	if (settings_.show_coverage_details && file_bam_!="")
	{
	//	writeCoverageReport(stream, file_bam_, file_roi_, roi_, genes_, settings_.min_depth, db_, settings_.recalculate_avg_depth, &roi_stats_, settings_.roi_low_cov);

		//writeCoverageReportCCDS(stream, file_bam_, genes_, settings_.min_depth, 0, db_, &roi_stats_, false, false);

		//writeCoverageReportCCDS(stream, file_bam_, genes_, settings_.min_depth, 5, db_, nullptr, true, true);
	}

	//OMIM table
	if (settings_.show_omim_table)
	{
		//prepare queries
		SqlQuery q_genes = db_.getQuery();
		q_genes.prepare("SELECT id, mim FROM omim_gene WHERE gene=:1");

		stream << "<p><b>" << trans("OMIM Gene und Phänotypen") << "</b>" << endl;
		stream << "</p>" << endl;
		stream << "<table>" << endl;
		stream << "<tr><td><b>" << trans("Gen") << "</b></td><td><b>" << trans("OMIM Gen MIM") << "</b></td><td><b>" << trans("OMIM Phänotypen") << "</b></td></tr>";
		foreach(QByteArray gene, genes_)
		{
			//approved gene symbol
			QByteArray gene_approved = db_.geneToApproved(gene, true);

			//generate table rows
			q_genes.bindValue(0, gene_approved);
			q_genes.exec();
			while (q_genes.next())
			{
				QString gene_id = q_genes.value(0).toByteArray();
				QString gene_mim = q_genes.value(1).toByteArray();

				stream << "<tr><td>" << gene << "</td><td>" << gene_mim << "</td><td>" << db_.getValues("SELECT phenotype FROM omim_phenotype WHERE omim_gene_id=" + gene_id).join("<br>")<< "</td></tr>";
			}
		}
		stream << "</table>" << endl;
	}

	//collect and display important tool versions
	if (settings_.show_tool_details)
	{
		stream << "<p><b>" << trans("Details zu Programmen der Analysepipeline") << "</b>" << endl;
		stream << "</p>" << endl;
		stream << "<table>" << endl;
		stream << "<tr><td><b>" << trans("Tool") << "</b></td><td><b>" << trans("Version") << "</b></td><td><b>" << trans("Parameter") << "</b></td></tr>";
		QStringList whitelist;
		whitelist << "SeqPurge" << "samblaster" << "/bwa" << "samtools" << "VcfLeftNormalize" <<  "freebayes" << "abra2" << "vcflib" << "ensembl-vep"; //current
		log_files_.sort();
		foreach(QString file, log_files_)
		{
			QStringList lines = Helper::loadTextFile(file);
			for(int i=0; i<lines.count(); ++i)
			{
				QString& line = lines[i];

				if (line.contains("Calling external tool")
					|| line.contains("command 1")
					|| line.contains("command 2")
					|| line.contains("command 3")
					|| line.contains("command 4")
					|| line.contains("command 5")
					|| line.contains("command 6")
					|| line.contains("command 7")
					|| line.contains("command 8")
					|| line.contains("command 9")
					|| line.contains("command 10")
					|| line.contains("command 11")
					|| line.contains("command 12")
					|| line.contains("command 13")
					|| line.contains("command 14")
					|| line.contains("command 15"))
				{
					//check if tool is whitelisted
					QString match = "";
					foreach(QString entry, whitelist)
					{
						if (line.contains(entry))
						{
							match = entry;
							break;
						}
					}
					if (match=="") continue;

					//extract version and parameters
					if (i+2>=lines.count()) continue;
					QString tool;
					if (line.contains("Calling external tool"))
					{
						tool = line.split("\t")[2].trimmed().mid(23, -1);
					}
					else
					{
						tool = line.split("\t")[2].trimmed().mid(13, -1);
					}
					tool.replace("/mnt/share/opt/", "[bin]/");
					QString version = lines[i+1].split("\t")[2].trimmed().mid(13);
					QString params = lines[i+2].split("\t")[2].trimmed().mid(13);
					params.replace(QRegExp("/tmp/[^ ]+"), "[file]");
					params.replace(QRegExp("\\./[^ ]+"), "[file]");
					params.replace("/mnt/share/data/", "[data]/");
					params.replace("//", "/");

					//output
					stream << "<tr><td>" << tool.toHtmlEscaped() << "</td><td>" << version.toHtmlEscaped() << "</td><td>" << params.toHtmlEscaped() << "</td></tr>" << endl;
				}
			}
		}
		stream << "</table>" << endl;
	}

	//close stream
	writeHtmlFooter(stream);
	outfile->close();


	validateAndCopyReport(temp_filename, file_rep_, true, false);

	//write XML file to transfer folder
	QString gsvar_variant_transfer = Settings::string("gsvar_variant_transfer");
	if (gsvar_variant_transfer!="")
	{
		QString xml_file = gsvar_variant_transfer + "/" + QFileInfo(file_rep_).fileName().replace(".html", ".xml");
		writeXML(xml_file, file_rep_);
	}
}

void ReportWorker::validateAndCopyReport(QString from, QString to, bool put_to_archive, bool is_rtf)
{
	//validate written file (HTML only)
	if (!is_rtf)
	{
		QString validation_error = XmlHelper::isValidXml(from);
		if (validation_error!="")
		{
			Log::warn("Generated HTML report at " + from + " is not well-formed: " + validation_error);
		}
	}

	if (QFile::exists(to) && !QFile(to).remove())
	{
		THROW(FileAccessException,"Could not remove previous " + QString(is_rtf ? "RTF" : "HTML") + " report: " + to);
	}
	if (!QFile::rename(from, to))
	{
		THROW(FileAccessException,"Could not move " + QString(is_rtf ? "RTF" : "HTML") + " report from temporary file " + from + " to " + to + " !");
	}

	//copy report to archive folder
	if(put_to_archive)
	{
		QString archive_folder = Settings::string("gsvar_report_archive");
		if (archive_folder!="")
		{
			QString file_rep_copy = archive_folder + "\\" + QFileInfo(to).fileName();
			if (QFile::exists(file_rep_copy) && !QFile::remove(file_rep_copy))
			{
				THROW(FileAccessException, "Could not remove previous " + QString(is_rtf ? "RTF" : "HTML") + " report in archive folder: " + file_rep_copy);
			}
			if (!QFile::copy(to, file_rep_copy))
			{
				THROW(FileAccessException, "Could not copy " + QString(is_rtf ? "RTF" : "HTML") + " report to archive folder: " + file_rep_copy);
			}
		}
	}
}

QByteArray ReportWorker::inheritance(const QByteArray& gene_info)
{
	QByteArrayList output;
	foreach(QByteArray gene, gene_info.split(','))
	{
		//extract inheritance info
		QByteArray inheritance;
		QByteArrayList parts = gene.replace('(',' ').replace(')',' ').split(' ');
		foreach(QByteArray part, parts)
		{
			if (part.startsWith("inh=")) inheritance = part.mid(4);
		}

		output << inheritance;
	}
	return output.join(",");
}

RtfSourceCode ReportWorker::trans(const QByteArray& text) const
{
	if (settings_.language=="german")
	{
		return text;
	}
	else if (settings_.language=="english")
	{
		QHash<QByteArray, QByteArray> de2en;
		de2en["Technischer Report zur bioinformatischen Analyse"] = "Technical Report for Bioinformatic Analysis";
		de2en["Probe"] = "Sample";
		de2en["Prozessierungssystem"] = "Processing system";
		de2en["Referenzgenom"] = "Reference genome";
		de2en["Datum"] = "Date";
		de2en["Benutzer"] = "User";
		de2en["Analysepipeline"] = "Analysis pipeline";
		de2en["Auswertungssoftware"] = "Analysis software";
		de2en["KASP-Ergebnis"] = " KASP result";
		de2en["Phänotyp"] = "Phenotype information";
		de2en["Filterkriterien"] = "Criteria for variant filtering";
		de2en["Gefundene Varianten in Zielregion gesamt"] = "Variants in target region";
		de2en["Anzahl Varianten ausgewählt für Report"] = "Variants selected for report";
		de2en["Varianten nach klinischer Interpretation im Kontext der Fragestellung"] = "List of prioritized variants";
		de2en["Vererbung"] = "Inheritance";
		de2en["Klasse"] = "Class";
		de2en["Details"] = "Details";
		de2en["Genotyp"] = "Genotype";
		de2en["Variante"] = "Variant";
		de2en["Gen"] = "Gene";
		de2en["Für Informationen zur Klassifizierung von Varianten, siehe allgemeine Zusatzinformationen."] = "For further information regarding the classification see Additional Information.";
		de2en["Teilweise können bei Varianten unklarer Signifikanz (Klasse 3) - in Abhängigkeit von der Art der genetischen Veränderung, der Familienanamnese und der Klinik des/der Patienten - weiterführende Untersuchungen eine änderung der Klassifizierung bewirken. Bei konkreten differentialdiagnostischen Hinweisen auf eine entsprechende Erkrankung ist eine humangenetische Mitbeurteilung erforderlich, zur Beurteilung ob erweiterte genetische Untersuchungen zielführend wären."] = "TODO";
		de2en["Klassifikation von Varianten"] = "Classification of variants";
		de2en["Die Klassifikation der Varianten erfolgt in Anlehnung an die Publikation von Plon et al. (Hum Mutat 2008)"] = "Classification and interpretation of variants: The classification of variants is based on the criteria of Plon et al. (PMID: 18951446). A short description of each class can be found in the following";
		de2en["Klasse 5: Eindeutig pathogene Veränderung / Mutation"] = "Class 5, pathogenic variant";
		de2en["Veränderung, die bereits in der Fachliteratur mit ausreichender Evidenz als krankheitsverursachend bezogen auf das vorliegende Krankheitsbild beschrieben wurde sowie als pathogen zu wertende Mutationstypen (i.d.R. Frameshift- bzw. Stoppmutationen)."] = "The variant is considered to be the cause of the patient's disease.";
		de2en["Klasse 4: Wahrscheinlich pathogene Veränderung"] = "Class 4, probably pathogenic variants";
		de2en["DNA-Veränderung, die aufgrund ihrer Eigenschaften als sehr wahrscheinlich krankheitsverursachend zu werten ist."] = "The identified variant is considered to be the probable cause of the patient's disease. This information should be used cautiously for clinical decision-making, as there is still a degree of uncertainty.";
		de2en["Klasse 3: Variante unklarer Signifikanz (VUS) - Unklare Pathogenität"] = "Class 3, variant of unclear significance (VUS)";
		de2en["Variante, bei der es unklar ist, ob eine krankheitsverursachende Wirkung besteht. Diese Varianten werden tabellarisch im technischen Report mitgeteilt."] = "The variant has characteristics of being an independent disease-causing mutation, but insufficient or conflicting evidence exists.";
		de2en["Klasse 2: Sehr wahrscheinlich benigne Veränderungen"] = "Class 2, most likely benign variants";
		de2en["Aufgrund der Häufigkeit in der Allgemeinbevölkerung oder der Lokalisation bzw. aufgrund von Angaben in der Literatur sehr wahrscheinlich benigne. Werden nicht mitgeteilt, können aber erfragt werden."] = "The variant is not likely to be the cause of the tested disease. Class 2 variants are not reported, but can be provided upon request.";
		de2en["Klasse 1: Benigne Veränderungen"] = "Class 1, benign variants";
		de2en["Werden nicht mitgeteilt, können aber erfragt werden."] = "The variant is not considered to be the cause of the tested disease. Class 1 variants are not reported, but can be provided upon request.";
		de2en["Zielregion"] = "Target region";
		de2en["Die Zielregion umfasst mindestens die CCDS (\"consensus coding sequence\") unten genannter Gene &plusmn;20 Basen flankierender intronischer Sequenz, kann aber auch zusätzliche Exons und/oder flankierende Basen beinhalten."] = "The target region includes CCDS (\"consensus coding sequence\") of the genes listed below &plusmn;20 flanking bases of the intronic sequence. It may comprise additional exons and/or flanking bases.";
		de2en["Name"] = "Name";
		de2en["Ausgewertete Gene"] = "Genes analyzed";
		de2en["OMIM Gene und Phänotypen"] = "OMIM gene and phenotypes";
		de2en["OMIM Phänotypen"] = "OMIM phenotypes";
		de2en["OMIM Gen MIM"] = "OMIM gene MIM";
		de2en["Gen"] = "Gene";
		de2en["Details zu Programmen der Analysepipeline"] = "Analysis pipeline tool details";
		de2en["Parameter"] = "Parameters";
		de2en["Version"] = "Version";
		de2en["Tool"] = "Tool";
		de2en["Abdeckungsstatistik"] = "Coverage statistics";
		de2en["Durchschnittliche Sequenziertiefe"] = "Average sequencing depth";
		de2en["Komplett abgedeckte Gene"] = "Genes without gaps";
		de2en["Anteil Regionen mit Tiefe <"] = "Percentage of regions with depth <";
		de2en["Fehlende Basen in nicht komplett abgedeckten Genen"] = "Number of missing bases for genes with gaps";
		de2en["Details Regionen mit Tiefe <"] = "Details regions with depth <";
		de2en["Koordinaten (hg19)"] = "Coordinates (hg19)";
		de2en["Chromosom"] = "Chromosome";
		de2en["Lücken"] = "Gaps";
		de2en["Abdeckungsstatistik für CCDS"] = "Coverage statistics for CCDS";
		de2en["Größe"] = "Size";
		de2en["Transcript"] = "Transcript";
		de2en["gesamt"] = "overall";
		de2en["mit Tiefe"] = "with depth";
		de2en["Geschlecht"] = "sample sex";
		de2en["Vater"] = "father";
		de2en["Mutter"] = "mother";


		if (!de2en.contains(text))
		{
			Log::warn("Could not translate to " + settings_.language + ": '" + text + "'");
		}

		return de2en[text];
	}

	THROW(ProgrammingException, "Unsupported language '" + settings_.language + "'!");
}

void ReportWorker::writeXML(QString outfile_name, QString report_document)
{
	QSharedPointer<QFile> outfile = Helper::openFileForWriting(outfile_name);

	QXmlStreamWriter w(outfile.data());
	w.setAutoFormatting(true);
	w.writeStartDocument();

	//element DiagnosticNgsReport
	w.writeStartElement("DiagnosticNgsReport");
	w.writeAttribute("version", "2");
	w.writeAttribute("type", settings_.report_type);

	//element ReportGeneration
	w.writeStartElement("ReportGeneration");
	w.writeAttribute("date", QDate::currentDate().toString("yyyy-MM-dd"));
	w.writeAttribute("user_name", Helper::userName());
	w.writeAttribute("software", QCoreApplication::applicationName() + " " + QCoreApplication::applicationVersion());
	w.writeAttribute("outcome", settings_.diag_status.outcome);
	w.writeEndElement();

	//element Sample
	w.writeStartElement("Sample");
	w.writeAttribute("name", sample_name_);

	SampleData sample_data = db_.getSampleData(db_.sampleId(sample_name_));
	w.writeAttribute("name_external", sample_data.name_external);
	ProcessedSampleData processed_sample_data = db_.getProcessedSampleData(db_.processedSampleId(sample_name_));
	w.writeAttribute("processing_system", processed_sample_data.processing_system);
	w.writeEndElement();

	//element TargetRegion (optional)
	if (file_roi_!="")
	{
		BedFile roi;
		roi.load(file_roi_);
		roi.merge();

		w.writeStartElement("TargetRegion");
		w.writeAttribute("name", QFileInfo(file_roi_).fileName().replace(".bed", ""));
		w.writeAttribute("regions", QString::number(roi.count()));
		w.writeAttribute("bases", QString::number(roi.baseCount()));
		QString gap_percentage = roi_stats_["gap_percentage"]; //cached from HTML report
		if (!gap_percentage.isEmpty())
		{
			w.writeAttribute("gap_percentage", gap_percentage);
		}
		QString ccds_sequenced = roi_stats_["ccds_sequenced"]; //cached from HTML report
		if (!ccds_sequenced.isEmpty())
		{
			w.writeAttribute("ccds_bases_sequenced", ccds_sequenced);
		}

		//contained genes
		foreach(const QString& gene, genes_)
		{
			w.writeStartElement("Gene");
			w.writeAttribute("name", gene);
			w.writeEndElement();
		}

		w.writeEndElement();
	}

	//element VariantList
	w.writeStartElement("VariantList");
	w.writeAttribute("overall_number", QString::number(variants_.count()));
	w.writeAttribute("genome_build", "hg19");

	//element Variant
	int geno_idx = variants_.getSampleHeader().infoByStatus(true).column_index;
	foreach(const ReportVariantConfiguration& var_conf, settings_.report_config.variantConfig())
	{
		if (var_conf.variant_type!=VariantType::SNVS_INDELS) continue;
		if (!var_conf.showInReport()) continue;
		if (var_conf.report_type!=settings_.report_type) continue;

		const Variant& variant = variants_[var_conf.variant_index];
		w.writeStartElement("Variant");
		w.writeAttribute("chr", variant.chr().str());
		w.writeAttribute("start", QString::number(variant.start()));
		w.writeAttribute("end", QString::number(variant.end()));
		w.writeAttribute("ref", variant.ref());
		w.writeAttribute("obs", variant.obs());
		w.writeAttribute("genotype", formatGenotype(processed_sample_data.gender.toLatin1(), variant.annotations()[geno_idx], variant));
		w.writeAttribute("causal", var_conf.causal ? "true" : "false");
		w.writeAttribute("de_novo", var_conf.de_novo ? "true" : "false");
		w.writeAttribute("comp_het", var_conf.comp_het ? "true" : "false");
		w.writeAttribute("mosaic", var_conf.mosaic ? "true" : "false");
		if (var_conf.inheritance!="n/a")
		{
			w.writeAttribute("inheritance", var_conf.inheritance);
		}
		QString classification = db_.getClassification(variant).classification;
		if (classification!="" && classification!="n/a")
		{
			w.writeAttribute("class", classification);
		}

		//element TranscriptInformation
		GeneSet genes;
		int i_co_sp = variants_.annotationIndexByName("coding_and_splicing", true, false);
		if (i_co_sp!=-1)
		{
			foreach(const VariantTranscript& trans, variant.transcriptAnnotations(i_co_sp))
			{
				w.writeStartElement("TranscriptInformation");
				w.writeAttribute("gene", trans.gene);
				w.writeAttribute("transcript_id", trans.id);
				w.writeAttribute("hgvs_c", trans.hgvs_c);
				w.writeAttribute("hgvs_p", trans.hgvs_p);
				w.writeAttribute("exon", trans.exon);
				w.writeEndElement();

				genes << trans.gene;
			}
		}

		//element Annotation
		for (int i=0; i<variant.annotations().count(); ++i)
		{
			if (i==geno_idx) continue;

			w.writeStartElement("Annotation");
			w.writeAttribute("name", variants_.annotations()[i].name());
			w.writeAttribute("value", variant.annotations()[i]);
			w.writeEndElement();
		}

		//element GeneDiseaseInformation
		if (var_conf.causal)
		{
			foreach(const QByteArray& gene, genes)
			{
				SqlQuery query = db_.getQuery();
				query.exec("SELECT dt.* FROM disease_gene dg, disease_term dt WHERE dt.id=dg.disease_term_id AND dg.gene='" + gene + "'");
				while(query.next())
				{
					w.writeStartElement("GeneDiseaseInformation");
					w.writeAttribute("gene", gene);
					w.writeAttribute("source", query.value("source").toString());
					w.writeAttribute("identifier", query.value("identifier").toString());
					w.writeAttribute("name", query.value("name").toString());
					w.writeEndElement();
				}
			}
		}

		//end of variant
		w.writeEndElement();
	}
	w.writeEndElement();

	//element ReportDocument
	w.writeStartElement("ReportDocument");
	QString format = QFileInfo(report_document).suffix().toUpper();
	w.writeAttribute("format", format);
	QByteArray base64_data = "";
	QFile file(report_document);
	file.open(QIODevice::ReadOnly);
	base64_data = file.readAll().toBase64();
	file.close();
	w.writeCharacters(base64_data);
	w.writeEndElement();

	w.writeEndDocument();
	outfile->close();

	//validate written XML file
	QString xml_error = XmlHelper::isValidXml(outfile_name, "://Resources/DiagnosticReport_v2.xsd");
	if (xml_error!="")
	{
		THROW(ProgrammingException, "ReportWorker::storeXML produced an invalid XML file: " + xml_error);
	}
}
