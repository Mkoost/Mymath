# Utilities

## Описание

Раздел посвящен пространству имен `mymath::utilities`. В нем содержаться различные функции упрощяющие работу с классами и со структурами Mymath.

## Краткие сведения

|Функция|описание|
|---|---|
| `print` | Выводит в `cout` содержание объекта |

##  Функции 

### `print(const &T)` 
Возможные реализации:
- `void print(const matrix&)`
- `void print(const vector&)`
- `void print(const quaternion&)`

Пример:
```c++
...
mymath::imat3 A{
    {1, 2, 3},
    {3, 1, 2},
    {2, 3, 1}};

mymath::utilities::print(A);
...
```

```terminal
1 2 3
3 1 2
2 3 1
```

### Cвязанное
- [Оглавление](./index.md)
- [Классы и структуры](./classes_and_structures.md)
- [Методы вычислений](./index.md)
