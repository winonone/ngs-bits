#ifndef FILTERDOCKWIDGET_H
#define FILTERDOCKWIDGET_H

#include <QDockWidget>
#include "ui_FilterDockWidget.h"
#include "VariantFilter.h"
#include "BedFile.h"

///Filter manager dock widget
class FilterDockWidget
	: public QDockWidget
{
	Q_OBJECT
	
public:
	/// Default constructor
	FilterDockWidget(QWidget *parent = 0);

	/// Resets to initial state (uncheck boxes, no ROI)
	void reset();
	/// Sets filter columns present in the open file
	void setFilterColumns(const QMap<QString, QString>& filter_cols);

    /// Applies predefined default filters (germline).
	void applyDefaultFilters();
    /// Applies predefined default filters (somatic).
    void applyDefaultFiltersSomatic();

	/// Returns if the MAF filter is enabled.
	bool applyMaf() const;
	/// Returns the maximum MAF percentage filter value.
	double mafPerc() const;

	/// Returns if the impact filter is enabled.
	bool applyImpact() const;
	/// Returns valid impact.
	QStringList impact() const;

	/// Returns if the classification filter is enabled.
	bool applyClassification() const;
	/// Returns the minimum classification.
	int classification() const;

	/// Returns if the genotype filter is enabled.
	bool applyGenotype() const;
	/// Returns the requested genotype.
	QString genotype() const;

	/// Returns if the IHDB filter is enabled.
	bool applyIhdb() const;
	/// Returns the maximum IHDB filter value.
	int ihdb() const;

	/// Returns if the trio filter is enabled.
	bool applyTrio() const;
	/// Returns if the compound-heterzygous filter is enabled.
	bool applyCompoundHet() const;

	/// Returns variants with a classification >= X should be kept (returns -1 if unset).
	int keepClassGreaterEqual() const;
	/// Returns variants with a classification 'M' (modifier) should be kept.
	bool keepClassM() const;

	///Returns the filter column terms to keep.
	QList<QByteArray> filterColumnsKeep() const;
	///Returns the filter column terms to remove.
	QList<QByteArray> filterColumnsRemove() const;
    ///Returns the filter column terms to filter.
    QList<QByteArray> filterColumnsFilter() const;

	/// Returns the target region file name or an empty string if unset.
	QString targetRegion() const;
	/// Returns the gene names filter.
	QList<QByteArray> genes() const;
	/// Returns the single target region filter.
	BedLine region() const;

	/// Returns the reference sample name or an empty string if unset.
	QString referenceSample() const;

	/// Returns a string representations of the applied filters
	QMap<QString, QString> appliedFilters() const;

signals:
	/// Signal that is emitted when an annotation filter changes its checkbox state, or if the ROI changes, or if the gene changes
	void filtersChanged();

protected slots:
	void addRoi();
	void addRoiTemp();
	void removeRoi();
	void roiSelectionChanged(int index);
	void addRef();
	void removeRef();
	void referenceSampleChanged(int index);
	void geneChanged();
	void regionChanged();
	void filterColumnStateChanged();

private:
	/// Loads the ROI filters from the INI file
	void loadROIFilters();
	/// Loads the reference file list of IGV
	void loadReferenceFiles();

    /// Resets the filters without blocking signals.
    void resetSignalsUnblocked();

	Ui::FilterDockWidget ui_;
	QList<QByteArray> last_genes_;
};

#endif // FILTERDOCKWIDGET_H
