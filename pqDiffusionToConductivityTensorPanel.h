// -*- c++ -*-

#ifndef __pqMaskFieldPanel_h
#define __pqMaskFieldPanel_h

// Qt Includes.
#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QWidget>
#include <QtDebug>

// ParaView Includes.
#include <pqAutoGeneratedObjectPanel.h>
#include <pqComboBoxDomain.h>
#include <vtkPVConfig.h>

class pqDiffusionToConductivityTensorPanel : public pqAutoGeneratedObjectPanel
{
  Q_OBJECT;
  typedef pqAutoGeneratedObjectPanel Superclass;

public:
  pqDiffusionToConductivityTensorPanel(pqProxy *proxy, QWidget *p);

private:


protected slots:

};

#endif //__pqMaskFieldPanel_h
