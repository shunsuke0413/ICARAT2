## 注意事項

リポジトリのルートディレクトリに存在する
CMakeLists.txt内の37行目を以下のように変換し，静的なビルドにしてください．

~~~
add_library(icarat STATIC ${SOURCES})
~~~
