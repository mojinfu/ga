package ga

import (
	"math/rand"
	"sort"
	"time"
)

type GeneAlgoInterface interface {
	// OrderMuta(*DNAIAndFitness)
	OrderMuta(*DNAIAndFitness)
	NodeMuta(int)
	Fitness(*DNAIAndFitness) float64
	Len() int
	// Get(int) interface{}
	// Key(int) string
	// Set(int, interface{})
	// Copy() GeneAlgoInterface
}

// type DNAI interface {
// 	SetOut(a interface{})
// 	GetDNAI(i int) int
// 	Len() int
// }

var iii int

func newDNAWithRandSort(this []int) []int {
	iii++
	out := make([]int, len(this))
	copy(out, this)
	for index := 0; index < len(out); index++ {
		newIndex := rand.Intn(len(out)-index) + index
		out[index], out[newIndex] = out[newIndex], out[index]
		// newIndexNode := newDna.Get(newIndex)
		// newDna.Set(newIndex, newDna.Get(index))
		// newDna.Set(index, newIndexNode)
		//	check(newDna)
	}
	return out
}
func SimpleOrderMuta(a *DNAIAndFitness, mutaP float64) {
	// newDnai := make([]int, dnai.Len())
	// for index := 0; index < dnai.Len(); index++ {
	// 	newDnai[index] = dnai.GetDNAI(index)
	// }
	var newIndex int
	for index := 0; index < len(a.dnai); index++ {
		if rand.Float64() < mutaP {
			newIndex = rand.Intn(len(a.dnai)-index) + index
			a.dnai[index], a.dnai[newIndex] = a.dnai[newIndex], a.dnai[index]
		}
	}
	a.fitness = nil
	a.out = nil
	return
}

// func SimpleOrderMuta(dnai []int, mutaP float64) []int {
// 	newDnai := make([]int, len(dnai))
// 	for index := 0; index < len(dnai); index++ {
// 		newDnai[index] = dnai[index]
// 	}
// 	for index := 0; index < len(dnai); index++ {
// 		if rand.Float64() < mutaP {
// 			newIndex := rand.Intn(len(newDnai)-index) + index
// 			newDnai[index], newDnai[newIndex] = newDnai[newIndex], newDnai[index]
// 		}
// 	}
// 	return newDnai
// }

type DNAIAndFitness struct {
	dnai    []int
	fitness *float64
	out     interface{}
}

func (this *DNAIAndFitness) GetOut() interface{} {
	return this.out
}
func (this *DNAIAndFitness) SetOut(a interface{}) {
	this.out = a
}
func (this *DNAIAndFitness) GetDNAI(i int) int {
	return this.dnai[i]
}
func (this *DNAIAndFitness) GetFitness(i int) float64 {
	if this.fitness == nil {
		return 0
	} else {
		return *this.fitness
	}
}
func (this *DNAIAndFitness) Len() int {
	return len(this.dnai)
}
func newDnaiAndFitness(a []int) *DNAIAndFitness {
	return &DNAIAndFitness{
		dnai:    a,
		fitness: nil,
	}
}

func (this *GA) Fitness(a *DNAIAndFitness) float64 {
	if a.fitness != nil {
		return *a.fitness
	} else {
		f := this.geneBank.Fitness(a)
		a.fitness = &f
		return *a.fitness
	}
}
func (this *GA) GA(bank GeneAlgoInterface) *DNAIAndFitness {
	this.initDNALength = bank.Len()
	this.geneBank = bank
	endAt := time.Now().Add(this.stopTrigger.timeLimit)
	initDnai := []int{}
	for i := 0; i < bank.Len(); i++ {
		initDnai = append(initDnai, i)
	}
	dnaList := []*DNAIAndFitness{
		{
			dnai:    initDnai,
			fitness: nil,
		},
	}

	var bestFitness float64 = bank.Fitness(dnaList[0])
	best := dnaList[0]
	for len(dnaList) < this.initialPopulationNum {
		dnaList = append(dnaList, newDnaiAndFitness(newDNAWithRandSort(best.dnai)))
	}
	// this.logGen(dnaList)
	for this.genNum = 1; this.genNum < this.stopTrigger.genNum; this.genNum++ {
		//cross
		if this.stopTrigger.timeLimit > 0 && endAt.Before(time.Now()) {
			break
		}
		afterCross := this.cross(dnaList)
		this.muta(afterCross)
		better := this.keepBetterDNA(dnaList)
		bestFitness = this.Fitness(better[0])
		// best = better[0]
		if this.Fitness(better[0]) > this.Fitness(best) {
			best = better[0]
		}
		// fmt.Println(this.genNum, "   this.fitness-->", this.Fitness(best), "------", this.Fitness(better[0]))
		if this.stopTrigger.bestFitnesstReach != nil && bestFitness >= *this.stopTrigger.bestFitnesstReach {
			break
		}
		newGen := []*DNAIAndFitness{}
		for len(newGen) < this.initialPopulationNum-len(better) {
			newGen = append(newGen, this.pickOne(afterCross))
		}
		// sed := (float64(this.genNum) / float64(this.stopTrigger.genNum))
		var sed int = 10
		if this.genNum%sed == 0 {
			if sed != 1 {
				sed--
			}
			newGen = append(newGen, better...)
		}
		if float64(this.genNum) > 0.2*float64(this.stopTrigger.genNum) {
			newGen = append(newGen, best)
		}
		dnaList = newGen
		// this.logGen(dnaList)
	}
	return best
}

type GA struct {
	initDNALength        int
	initialPopulationNum int
	genealogy            genealogyStruct
	pepoleMutaP          float64
	pepoleNodeMutaP      float64
	recordGenealogy      bool
	stopTrigger          stopTrigger
	Counter              Counter
	geneBank             GeneAlgoInterface

	genNum int
}
type genealogyStruct struct {
	dnaG            [][]GeneAlgoInterface
	bestFitnessList []float64
}
type Counter struct {
	genNum int
	Muta   int
}
type stopTrigger struct {
	timeLimit         time.Duration
	genNum            int
	bestFitnesstReach *float64
}

func NewGA() *GA {
	return &GA{
		// initDNA
		initDNALength:        -1,
		initialPopulationNum: 30,
		pepoleMutaP:          0.2,
		pepoleNodeMutaP:      0.2,
		stopTrigger: stopTrigger{
			timeLimit:         time.Millisecond * 1000 * 10,
			genNum:            30,
			bestFitnesstReach: nil,
		},
		genNum: 0,
	}
}

func (this *GA) SetPepoleMutaP(a float64) {
	this.pepoleMutaP = a
}
func (this *GA) FitnessGenealogy() []float64 {
	return this.genealogy.bestFitnessList
}
func (this *GA) SetBestFitnessNeed(f float64) {
	this.stopTrigger.bestFitnesstReach = &f
}
func (this *GA) SetTimeLimit(l time.Duration) {
	this.stopTrigger.timeLimit = l
}
func (this *GA) SetMaxGenerNum(m int) {
	this.stopTrigger.genNum = m
}
func (this *GA) IfRecordGenealogy(r bool) {
	this.recordGenealogy = r
}
func (this *GA) GetGenNum() int {
	return this.genNum
}
func (this *GA) GetGenNumLimit() int {
	return this.stopTrigger.genNum
}

// func (this *GA) logGen(dnaList []*DNAIAndFitness) {
// 	if len(dnaList) < 0 {
// 		panic("")
// 	}
// 	this.Counter.genNum++
// 	if this.recordGenealogy {
// 		this.genealogy.dnaG = append(this.genealogy.dnaG, dnaList)
// 	}
// 	this.genealogy.bestFitnessList = append(this.genealogy.bestFitnessList, dnaList[0].Fitness())
// }

type sortDNAList struct {
	list []*DNAIAndFitness
	ga   *GA
}

func (this sortDNAList) Less(a, b int) bool {
	return this.ga.Fitness(this.list[a]) > this.ga.Fitness(this.list[b])
}
func (this sortDNAList) Swap(a, b int) {
	this.list[a], this.list[b] = this.list[b], this.list[a]
}
func (this sortDNAList) Len() int {
	return len(this.list)
}
func (this *GA) keepBetterDNA(dnaList []*DNAIAndFitness) []*DNAIAndFitness {
	sort.Sort(sortDNAList{
		list: dnaList,
		ga:   this,
	})
	var newGen []*DNAIAndFitness
	for i := 0; i < 1+len(dnaList)/10 && i < len(dnaList); i++ {
		newGen = append(newGen, dnaList[i])
	}
	return newGen
}
func (this *GA) cross(dnaList []*DNAIAndFitness) []*DNAIAndFitness {
	sort.Sort(sortDNAList{
		list: dnaList,
		ga:   this,
	})
	var newGen []*DNAIAndFitness
	for {
		if len(newGen) < int(float64(len(dnaList))*1.5) {
			t := this.pickTwoDifferent_old(dnaList)
			list := crossAAndB(t[0].dnai, t[1].dnai)
			for i := range list {
				newGen = append(newGen, newDnaiAndFitness(list[i]))
			}

		} else {
			break
		}
	}
	return newGen
}

var corssRank int

func crossAAndB(a, b []int) [][]int {
	// a.(*baseInfoDNAi).CC()
	// b.(*baseInfoDNAi).CC()
	corssRank++
	var c []int = make([]int, len(a))
	unionMap := make(map[int]bool)
	bitset := make([]bool, len(a))

	for i := 0; i < len(a); i++ {
		if rand.Float64() > 0.5 {
			c[i] = a[i] //c.Set(i, a.Get(i)) //
			unionMap[a[i]] = true
			bitset[i] = true
		}
	}
	bi := 0
	for ci := 0; ci < len(a); ci++ {
		if bitset[ci] {
			continue
		}
		for bi < len(b) {
			_, isOK := unionMap[b[bi]]
			if !isOK {
				unionMap[b[bi]] = true
				c[ci] = b[bi]
				bi++
				break
			} else {
				bi++
			}
		}
	}
	//for i := range newGen {
	// fmt.Println("corssRank:", corssRank)
	//	c.(*baseInfoDNAi).CC()
	//}
	//	check(c)
	return [][]int{c}

}
func (this *GA) muta(dnaiList []*DNAIAndFitness) {
	sort.Sort(sortDNAList{
		list: dnaiList,
		ga:   this,
	})
	for i := 0; i < len(dnaiList); i++ {
		if i < 1+len(dnaiList)/5 {
			continue
		}
		if rand.Float64() < this.pepoleMutaP {
			this.Counter.Muta++
			this.geneBank.OrderMuta(dnaiList[i]) //todo
			dnaiList[i].fitness = nil
			dnaiList[i].out = nil
		}
	}
}

func (this *GA) pickOne(dnaList []*DNAIAndFitness) *DNAIAndFitness {
	a := dnaList[rand.Intn(len(dnaList))]
	b := dnaList[rand.Intn(len(dnaList))]

	if this.Fitness(a) > this.Fitness(b) {
		return a
	} else {
		return b
	}
}

func (this *GA) pickTwoDifferent(dnaList []*DNAIAndFitness) [2]*DNAIAndFitness {
	if len(dnaList) < 2 {
		panic("")
	}
	if len(dnaList) == 2 {
		return [2]*DNAIAndFitness{dnaList[0], dnaList[1]}
	}
	l := len(dnaList)
	var a, b, c int
	a = rand.Intn(l)
	bI := rand.Intn(l-2) + 1
	b = (a + bI) % l
	abs := a - b
	if abs < 0 {
		abs = -abs
	}
	if a-b == 1 || b-a == 1 {
		return [2]*DNAIAndFitness{dnaList[a], dnaList[b]}
	}
	cI := rand.Intn(abs-1) + 1
	if a < b {
		c = a + cI
	} else {
		c = b + cI
	}

	// for {
	// 	c = rand.Intn(l)
	// 	if c != a && c != b {
	// 		break
	// 	}
	// }
	temp := []*DNAIAndFitness{dnaList[a], dnaList[b], dnaList[c]}
	sort.Sort(sortDNAList{
		list: temp,
		ga:   this,
	})
	return [2]*DNAIAndFitness{temp[0], temp[1]}

	//return [2]*DNAIAndFitness{dnaList[a], dnaList[b]}
}

func (this *GA) pickTwoDifferent_old(dnaList []*DNAIAndFitness) [2]*DNAIAndFitness {
	if len(dnaList) < 2 {
		panic("")
	}
	if len(dnaList) == 2 {
		return [2]*DNAIAndFitness{dnaList[0], dnaList[1]}
	}
	var a, b, c int
	a = rand.Intn(len(dnaList))

	for {
		b = rand.Intn(len(dnaList))
		if b != a {
			break
		}
	}
	for {
		c = rand.Intn(len(dnaList))
		if c != a && c != b {
			break
		}
	}
	temp := []*DNAIAndFitness{dnaList[a], dnaList[b], dnaList[c]}
	sort.Sort(sortDNAList{
		list: temp,
		ga:   this,
	})
	return [2]*DNAIAndFitness{temp[0], temp[1]}
}
