FetchContent_Declare(
  efika-core
  GIT_REPOSITORY git@github.com:jiverson002/efika-core.git
)

FetchContent_Declare(
  efika-data
  GIT_REPOSITORY git@github.com:jiverson002/efika-data.git
)

FetchContent_Declare(
  efika-io
  GIT_REPOSITORY git@github.com:jiverson002/efika-io.git
)

FetchContent_MakeAvailable(efika-core efika-data efika-io)
