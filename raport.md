# Analiza rakowego środowiska immunologicznego

## Opis

Celem projektu jest stworzenie programu do przenalizowania danych z próbek pacjentów, między innymi w poszukiwaniu TLSów. TLS (Tertiary Lymphoid Structures) to struktury tkanki chłonnej. Potencjalnym sygnałem, że w jakimś miejscu znajduje się TLS, może być znaczna ilość komórek B w danym miejscu. Obecność TLSów ma istotne znaczenie kliniczne, zwłaszcza w kontekście odpowiedzi immunologicznej i prognozy choroby.  

## Parametry

Użytkownik ma do wyboru 4 interaktywnie wybierane parametry:

- Wybór pliku, z którego ciągniemy dane
- Radius - na podstawie którego tworzony jest graf sąsiedztwa Bcell'i. 
- Minimalna ilość Bcell'i, żeby komponent był brany pod uwagę w dalszej analizie
- Component który chcemy w danym momencie zwizualizować

## Opis metody

Tworzę grafy sąsiedztwa komórek B z zadanym promieniem. Później dzielę te dane na spójne componenty, a wreszcie filtruje je zgodnie z zadaną przez użytkownika ilością komórek B, które muszą znajdować się w danym componencie aby był brany pod uwagę w analizie.
*Uwaga* - Jeśli program nie znajdzie żadnego componentu z conajmniej zadaną liczbą Bcell'i, to zwróci component z największą liczbą Bcell'i jaki znalazł.

Następnie, dla wszystkich komórek B znajdujących się w danym componencie, szukam sąsiadów, przez co ostatecznie dostaję sąsiedztwo całego componentu (oczywiście usuwając zdublowanych sąsiadów).

Z otrzymanych danych dotyczących sąsiedztwa componentów dla każdego z nich liczę procentowy udział typów komórek w ich otoczeniu. Na podstawie tych danych liczę odległości między każdą parą componentów, i tworzę dendrogram użwając do tego klastrowania hierarchicznego.

Później użytkownik może wybrać, informację o którym komponencie w danej chwili chciałby zobaczyć.
Informację te obejmują ilość komórek B w danym componencie, oraz histogram z procentowym udziałem typów komórek. Na podstawie dendrogramu użytkownik jest w stanie zobaczyć, które z danych componentów są do siebie najbardziej (lub najmniej) podobne (jeśli chodzi o procentowy udział typów komórek w otoczeniu danego componentu). Dane te mogą być wykorzystywane do szukania kandydatów na TLS. 
