# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PolygonByPolarCoordinates
                                 A QGIS plugin
 Creates polygon deifiedn by vertices in polar coordinates system
                             -------------------
        begin                : 2018-08-05
        copyright            : (C) 2018 by Paweł Strzleeiwcz
        email                : @
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load PolygonByPolarCoordinates class from file PolygonByPolarCoordinates.

    :param iface: A QGIS interface instance.
    :type iface: QgisInterface
    """
    #
    from .polygon_by_polar_coordinates import PolygonByPolarCoordinates
    return PolygonByPolarCoordinates(iface)
