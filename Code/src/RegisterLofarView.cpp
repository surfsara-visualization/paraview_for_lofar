/*=========================================================================

   Program: ParaView
   Module:    RegisterLofarView.h

   Copyright (c) 2005,2006 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2.

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

========================================================================*/
#include "RegisterLofarView.h"

#include "pqApplicationCore.h"
#include "pqInterfaceTracker.h"
#include "pqStandardViewModules.h"
#include "vtkSMSessionProxyManager.h"
#include "pqServer.h"

#include <QStringList>

namespace
{
  /// interface class for plugins that create view modules
  class VTK_EXPORT LofarViewInterface : public pqStandardViewModules
  {
public:
  LofarViewInterface(QObject* o) : pqStandardViewModules(o) { }
  ~LofarViewInterface() { }

  virtual QStringList viewTypes() const
    {
    QStringList list;
    list << "LofarView";
    return list;
    }

  QString viewTypeName(const QString& val) const
    {
    if (val == "LofarView")
      {
      return "LofarView";
      }
    return "LofarView";
    }

  virtual vtkSMProxy* createViewProxy(
    const QString& viewtype, pqServer *server)
    {
    vtkSMSessionProxyManager* pxm = server->proxyManager();
    if (viewtype == "LofarView")
      {
      return pxm->NewProxy("views", "LofarView");
      }
    return NULL;
    }
  };
}

//-----------------------------------------------------------------------------
RegisterLofarView::RegisterLofarView(QObject* p)
  : Superclass(p)
{
}

//-----------------------------------------------------------------------------
RegisterLofarView::~RegisterLofarView()
{
}


//-----------------------------------------------------------------------------
void RegisterLofarView::startup()
{
  // Register ParaView interfaces.
  pqInterfaceTracker* pgm = pqApplicationCore::instance()->interfaceTracker();

  // * adds support for standard paraview views.
  pgm->addInterface(new LofarViewInterface(pgm));
}


//-----------------------------------------------------------------------------
void RegisterLofarView::shutdown()
{
}
