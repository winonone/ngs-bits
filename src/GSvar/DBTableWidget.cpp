#include "DBTableWidget.h"
#include "Exceptions.h"
#include "GUIHelper.h"

#include <QHeaderView>
#include <QAction>
#include <QApplication>
#include <QClipboard>
#include <QKeyEvent>

DBTableWidget::DBTableWidget(QWidget* parent)
	: QTableWidget(parent)
{
	//general settings
	verticalHeader()->setVisible(false);
	setSelectionMode(QAbstractItemView::ExtendedSelection);
	setSelectionBehavior(QAbstractItemView::SelectRows);
	setWordWrap(false);

	//context menu
	setContextMenuPolicy(Qt::ActionsContextMenu);
	QAction* copy_action = new QAction(QIcon(":/Icons/CopyClipboard.png"), "Copy", this);
	addAction(copy_action);
	connect(copy_action, SIGNAL(triggered(bool)), this, SLOT(copyToClipboard()));
}

void DBTableWidget::setData(const DBTable& table)
{
	QStringList headers = table.headers();

	//table
	table_ = table.tableName();

	//resize
	clearContents();
	setRowCount(table.rowCount());
	setColumnCount(headers.count());

	//headers
	for(int c=0; c<headers.count(); ++c)
	{
		setHorizontalHeaderItem(c, createItem(headers[c], Qt::AlignCenter));
	}

	//content
	ids_.clear();
	for(int r=0; r<table.rowCount(); ++r)
	{
		const DBRow& row = table.row(r);
		ids_ <<row.id();
		for(int c=0; c<headers.count(); ++c)
		{
			setItem(r, c, createItem(row.value(c)));
		}
	}

	//fomatting
	GUIHelper::resizeTableCells(this, 200);
}

int DBTableWidget::columnIndex(const QString& column_header) const
{
	for (int c=0; c<columnCount(); ++c)
	{
		if (horizontalHeaderItem(c)->text()==column_header)
		{
			return c;
		}
	}

	THROW(ArgumentException, "Could not find column with header '" + column_header + "'");
}

void DBTableWidget::setColumnTooltips(const QString& column_header, const QStringList& tooltips)
{
	if (tooltips.count()!=rowCount())
	{
		THROW(ArgumentException, "Invalid tooltip count '" + QString::number(tooltips.count()) + "' in DBTableWidget::setColumnTooltips - expected '" + QString::number(rowCount()) + "'!");
	}

	int c = columnIndex(column_header);
	for(int r=0; r<rowCount(); ++r)
	{
		QTableWidgetItem* table_item = item(r, c);
		if (table_item==nullptr) continue;
		table_item->setToolTip(tooltips[r]);
	}
}

void DBTableWidget::setColumnColors(const QString& column_header, const QList<QColor>& colors)
{
	if (colors.count()!=rowCount())
	{
		THROW(ArgumentException, "Invalid color count '" + QString::number(colors.count()) + "' in DBTableWidget::setColumnColors - expected '" + QString::number(rowCount()) + "'!");
	}

	int c = columnIndex(column_header);
	for (int r=0; r<rowCount(); ++r)
	{
		const QColor& color = colors[r];
		if (!color.isValid()) continue;

		QTableWidgetItem* table_item = item(r, c);
		if (table_item==nullptr) continue;

		table_item->setBackgroundColor(color);
	}
}

void DBTableWidget::setBackgroundColorIfEqual(const QString& column_header, const QColor& color, const QString& text)
{
	setBackgroundColorIf(column_header, color, [text](const QString& str) { return str==text; });
}

void DBTableWidget::setBackgroundColorIfLt(const QString& column_header, const QColor& color, double cutoff)
{
	setBackgroundColorIf(column_header, color, [cutoff](const QString& str) { bool ok; double value = str.toDouble(&ok); if (!ok) return false; return value<cutoff; });
}

void DBTableWidget::setBackgroundColorIfGt(const QString& column_header, const QColor& color, double cutoff)
{
	setBackgroundColorIf(column_header, color, [cutoff](const QString& str) { bool ok; double value = str.toDouble(&ok); if (!ok) return false; return value>cutoff; });
}

QSet<int> DBTableWidget::selectedRows() const
{
	QSet<int> output;

	foreach(const QTableWidgetSelectionRange& range, selectedRanges())
	{
		for (int row=range.topRow(); row<=range.bottomRow(); ++row)
		{
			output << row;
		}
	}

	return output;
}

const QString& DBTableWidget::getId(int r) const
{
	if (r<0 || r>=rowCount())
	{
		THROW(ArgumentException, "Invalid row index '" + QString::number(r) + "' in DBTableWidget::getId!");
	}

	return ids_[r];
}

const QString& DBTableWidget::tableName() const
{
	return table_;
}

QTableWidgetItem* DBTableWidget::createItem(const QString& text, int alignment)
{
	auto item = new QTableWidgetItem();
	item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);

	item->setText(text);

	item->setTextAlignment(alignment);

	return item;
}

void DBTableWidget::keyPressEvent(QKeyEvent* event)
{
	if(event->matches(QKeySequence::Copy))
	{
		copyToClipboard();
		event->accept();
		return;
	}

	QTableWidget::keyPressEvent(event);
}

void DBTableWidget::copyToClipboard()
{
	//header
	QString output = "#";
	for (int col=0; col<columnCount(); ++col)
	{
		if (col!=0) output += "\t";
		output += horizontalHeaderItem(col)->text();
	}
	output += "\n";

	//rows
	QSet<int> selected_rows = selectedRows();
	for (int row=0; row<rowCount(); ++row)
	{
		//skip hidden
		if (isRowHidden(row)) continue;

		//skip unselected
		if (selected_rows.count()>0 && !selected_rows.contains(row)) continue;

		for (int col=0; col<columnCount(); ++col)
		{
			if (col!=0) output += "\t";
			output += item(row, col)->text().replace('\t', ' ').replace('\n', ' ');
		}
		output += "\n";
	}

	QApplication::clipboard()->setText(output);
}
