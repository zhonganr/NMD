/****************************************************************************
** Meta object code from reading C++ file 'GUI.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../include/GUI.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GUI.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_GUI_t {
    QByteArrayData data[13];
    char stringdata0[269];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_GUI_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_GUI_t qt_meta_stringdata_GUI = {
    {
QT_MOC_LITERAL(0, 0, 3), // "GUI"
QT_MOC_LITERAL(1, 4, 21), // "SetPage_CalculateMode"
QT_MOC_LITERAL(2, 26, 0), // ""
QT_MOC_LITERAL(3, 27, 24), // "Index_Page_CalculateMode"
QT_MOC_LITERAL(4, 52, 20), // "SetPage_ElasticField"
QT_MOC_LITERAL(5, 73, 23), // "Index_Page_ElasticField"
QT_MOC_LITERAL(6, 97, 17), // "Pthread_State_Run"
QT_MOC_LITERAL(7, 115, 30), // "Pthread_State_Run_Part_Prepare"
QT_MOC_LITERAL(8, 146, 26), // "Pthread_State_Run_Part_Run"
QT_MOC_LITERAL(9, 173, 33), // "Pthread_State_Run_Part_Run_St..."
QT_MOC_LITERAL(10, 207, 18), // "Pthread_State_Stop"
QT_MOC_LITERAL(11, 226, 22), // "Pthread_State_Continue"
QT_MOC_LITERAL(12, 249, 19) // "Pthread_State_Break"

    },
    "GUI\0SetPage_CalculateMode\0\0"
    "Index_Page_CalculateMode\0SetPage_ElasticField\0"
    "Index_Page_ElasticField\0Pthread_State_Run\0"
    "Pthread_State_Run_Part_Prepare\0"
    "Pthread_State_Run_Part_Run\0"
    "Pthread_State_Run_Part_Run_Strain\0"
    "Pthread_State_Stop\0Pthread_State_Continue\0"
    "Pthread_State_Break"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_GUI[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   59,    2, 0x08 /* Private */,
       4,    1,   62,    2, 0x08 /* Private */,
       6,    0,   65,    2, 0x08 /* Private */,
       7,    0,   66,    2, 0x08 /* Private */,
       8,    0,   67,    2, 0x08 /* Private */,
       9,    0,   68,    2, 0x08 /* Private */,
      10,    0,   69,    2, 0x08 /* Private */,
      11,    0,   70,    2, 0x08 /* Private */,
      12,    0,   71,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void GUI::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<GUI *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->SetPage_CalculateMode((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->SetPage_ElasticField((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->Pthread_State_Run(); break;
        case 3: _t->Pthread_State_Run_Part_Prepare(); break;
        case 4: _t->Pthread_State_Run_Part_Run(); break;
        case 5: _t->Pthread_State_Run_Part_Run_Strain(); break;
        case 6: _t->Pthread_State_Stop(); break;
        case 7: _t->Pthread_State_Continue(); break;
        case 8: _t->Pthread_State_Break(); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject GUI::staticMetaObject = { {
    QMetaObject::SuperData::link<QMainWindow::staticMetaObject>(),
    qt_meta_stringdata_GUI.data,
    qt_meta_data_GUI,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *GUI::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *GUI::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_GUI.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int GUI::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 9)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 9;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
