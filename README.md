# sagd_new
### Что тут вообще
Решается задача одномерной трехфазной фильтрации.
[Уравнения](https://drive.google.com/file/d/0B7GsTw5yAyJbb1F0OVJEbU5FWEk/view?usp=sharing)
Пока есть три этапа:

1. Трехфазная с противопоточной аппроксимацией (Трехфазная вынесена в отдельную ветку)
2. Двухфазная (скелет неподвижен, постоянная пористость) с решением задачи о распаде разрыва.
3. Двухфазная (скелет неподвижен, постоянная пористость) с решением упрощенной задачи о распаде разрыва - HLL. Могут быть
только 2 волны, без контактных разрывов

---------------------------------------------------------------------------------------------------------------------

Плюсом еще выяснял гиперболическая ли задача (та что с тремя фазами), получилось что вообще то не везде гиперболическая.
